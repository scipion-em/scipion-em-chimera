# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Roberto Marabini 
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************


from os.path import exists
import os
import time
import json
import pyworkflow.protocol.params as params
import requests
import xml.etree.ElementTree as ET
import glob
from pwem.objects.data import Sequence
import pwem.convert as emconv
from pwem.convert import AtomicStructHandler
import pwem.objects as emobj

from pwem.protocols import EMProtocol
from pwem.viewers.viewer_chimera import chimeraScriptFileName, Chimera
from chimera import Plugin
from chimera.utils import getEnvDictionary
from pwem.convert import AtomicStructHandler
from chimera.colabs.browser import createColabScript

class ProtImportAtomStructAlphafold(EMProtocol):
    """ Protocol to import atomic structures generated by alphafold.\n
    If you choose the "Execute alphafold Locally" option you will need 
    a local alphafold NO docker instalation as
    described here: https://github.com/kalininalab/alphafold_non_docker
    """
    _label = 'alphafold prediction'
    # SEQUENCEFILENAME = '_sequence.fasta'
    IMPORT_FROM_EBI = 0
    IMPORT_FROM_SEQ_BLAST = 1
    IMPORT_REMOTE_ALPHAFOLD = 2
    IMPORT_LOCAL_ALPHAFOLD = 3
    INPUTFASTAFILE = 'seqs'    

    CHIMERA = 0
    PHENIX = 1
    TEST = 2

    url = {}
    url[CHIMERA] = "https://colab.research.google.com/github/scipion-em/scipion-em-chimera/blob/devel/chimera/colabs/chimera_alphafold_colab.ipynb"
    url[PHENIX]  = "https://colab.research.google.com/github/scipion-em/scipion-em-chimera/blob/devel/chimera/colabs/phenix_alphafold_colab.ipynb"
    url[TEST]  = "https://colab.research.google.com/github/scipion-em/scipion-em-chimera/blob/devel/chimera/colabs/test_colab.ipynb"


    def __init__(self, **args):
        EMProtocol.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('source', params.EnumParam,
                      choices=['Chimera', 
                               'Phenix', 
                               'colabfold', 
                               'local alphafold',
                               ],
                      display=params.EnumParam.DISPLAY_HLIST,
                      label="Source ",
                      default=self.IMPORT_FROM_EBI,
                      help='Search alphafold model in:\n '
                           '* EBI database by uniprot ID\n '
                           '* EBI database by homologous (Blast)\n '
                           '* Execute alphafold in Google-colab'
                           '* Execute alphafold Locally (multimer supported)\n')
        form.addParam('uniProtId', params.StringParam,
                      condition='source == %d'  %
                                (self.IMPORT_FROM_EBI),
                      label="UniProt name/ID ", allowsNull=True,
                      help='Write a UniProt ID (six or ten alphanumeric '
                           'characters; examples: A2BC19, P12345, '
                           'A0A022YWF9, DGAL_ECOLI).\n You can convert other '
                           'database identifiers to UniProt accession codes '
                           'by using the "ID Mapping" tab on '
                           'https://www.uniprot.org/')
        # TODO: add setofsequences to colab
        # so far this version only handles a single sequence
        form.addParam('inputSequence', params.PointerParam, pointerClass="Sequence",
                       label='Reference sequence', allowsNull=True,
                       condition='source == %d or '
                                 'source == %d '  % (self.IMPORT_FROM_SEQ_BLAST,
                                                     self.IMPORT_REMOTE_ALPHAFOLD),
                       help="Input the aminoacid sequence to blast or send to colab lab")
        # # list different colabs if source == IMPORT_REMOTE_ALPHAFOLD
        form.addParam('colabID', params.EnumParam,
                        choices=['Chimera', 
                                 'Phenix', 
                                 'Test'
                                 ],
                        display=params.EnumParam.DISPLAY_HLIST,
                        label="Colab Notebook ",
                        default=self.CHIMERA,
                        condition='source == %d ' % self.IMPORT_REMOTE_ALPHAFOLD,
                        help='Execute alphafold in Google-colab.\n'
                            '  Two notebooks are available from\n'
                             'Chimera and Phenix respectively'
                             )
        form.addParam('useTemplatesFromPDB', params.IntParam,
                      label='Use Templates from PDB',
                      default=-1,
                      condition='source == %d and colabID == %d' % (self.IMPORT_REMOTE_ALPHAFOLD, self.PHENIX),
                      help="Number of PDB templates to use. Set to -1 to disable"
                           "Suggested value=20"
                    )
                    
        form.addParam('inputSequenceS', params.MultiPointerParam,
                      pointerClass="Sequence", allowsNull=True,
                      label='Structures',
                      condition='source == %d '  % (self.IMPORT_LOCAL_ALPHAFOLD),
                      help="Structures to be procesed by local Alphafold ")

        form.addParam('maxTemplateDate', params.StringParam,
                      label='Use Template until',
                      default='2050-01-01',
                      condition='source == %d' % (self.IMPORT_LOCAL_ALPHAFOLD),
                      help="Maximum template release date to consider (YYYY-MM-DD)")

        form.addParam('isProkaryote', params.BooleanParam,
                      label='Is prokaryote?',
                      default=False,
                      condition='source == %d' % (self.IMPORT_LOCAL_ALPHAFOLD),
                      help="Optional for multimer system, not used by the single "
                           "chain system. A boolean specifying true where the target "
                           "complex is from a prokaryote, and false where it is not, "
                           "or where the origin is unknown. This value determine the "
                           "pairing method for the MSA (default: 'None')"
                        )
        form.addParam('doGpu', params.BooleanParam, default=True,
                      condition='source == %d' % (self.IMPORT_LOCAL_ALPHAFOLD),
                      label='Use GPU acceleration?',
                      help='If set to Yes, the job will try to use GPU '
                           'acceleration.')
        form.addParam('gpusToUse', params.StringParam,
                      label='Which GPUs to use:',
                      condition='source == %d and doGpu' % (self.IMPORT_LOCAL_ALPHAFOLD),
                      default = '0',
                      help='This argument is not necessary. If left empty, '
                           'the job itself will try to allocate available '
                           'GPU resources. You can override the default '
                           'allocation by providing a list of which GPUs '
                           '(0,1,2,3, etc) to use.')

        form.addParam('extraFlags', params.StringParam,
                          default = '',
                          condition='source == %d' % (self.IMPORT_LOCAL_ALPHAFOLD),
                          label = 'Extra flags',
                          help = "# -n <openmm_threads>   OpenMM threads (default: all available cores) "
                                 "-c <db_preset>        Choose preset MSA database configuration - smaller "
                                 "genetic database config (reduced_dbs) or full genetic database config "
                                 "(full_dbs) (default: 'full_dbs')"
                    )

        form.addParam('extraCommands', params.StringParam,
                          default = '',
                          condition = 'False',
                          label = 'Extra commands for alphafold',
                          help = "Add extra commands in cmd file. Use for testing")
        form.addParam('hideMessage', params.BooleanParam, default=False,
                      condition='source != %d and source!= %d' % (self.IMPORT_FROM_EBI, self.IMPORT_LOCAL_ALPHAFOLD),
                      label='Hide help popup window',
                      help='If set to Yes no help message will be shown in chimera at start up.')

    def _getDefaultParallel(self):
        """This protocol doesn't have mpi version"""
        return (0, 0)


    def _insertAllSteps(self):
        hideMessage = self.hideMessage.get()
        if self.source == self.IMPORT_FROM_EBI:
            uniProtId = self.uniProtId.get()
            self._insertFunctionStep('_getModelFromEBI', uniProtId)
        elif self.source == self.IMPORT_FROM_SEQ_BLAST:
            inputSequence = self.inputSequence.get().getSequence()
            self._insertFunctionStep('_getModelFromBlast', inputSequence, hideMessage)
        elif self.source == self.IMPORT_REMOTE_ALPHAFOLD:
            inputSequence = self.inputSequence.get().getSequence()
            colabID = self.colabID.get()
            useTemplatesFromPDB = self.useTemplatesFromPDB.get()
            self._insertFunctionStep('_getModelFromColab', inputSequence, colabID, hideMessage, useTemplatesFromPDB)
        elif self.source == self.IMPORT_LOCAL_ALPHAFOLD:
            seqs = []
            for seq in self.inputSequenceS:
                s = seq.get()
                seqs.append( (s.getId(), s.getSequence() ) )
            isProkaryote = self.isProkaryote.get()
            self._insertFunctionStep('_getModelFromLocal', seqs, isProkaryote, 
                                     self.maxTemplateDate.get(),
                                     self.doGpu.get(), self.gpusToUse.get(),
                                     self.extraFlags.get())
        else:
            print("WRONG source")

    def createInputFastaFile(self, seqs):
        """Get sequence as string and create the corresponding fata file"""
        fastaFileName = self._getExtraPath(self.INPUTFASTAFILE + ".fasta")
        f = open(fastaFileName, "w")
        for id, seq in seqs:
            f.write(f"> {id}\n")
            f.write(f"{seq}\n")
        f.close()
        return os.path.abspath(fastaFileName)


    def _getModelFromLocal(self, seqs, isProkaryote, maxTemplateDate, doGpu, gpusToUse, extraFlags):
        """Use local alphafold installed in a conda enviroment
        variable CONDA_ACTIVATION_CMD and ALPHAFOLD_HOME are needed
        """
        # get environ variables.
        CONDA_ACTIVATION_CMD = Plugin.getCondaActivationCmd()
        CONDA_ACTIVATION_CMD = CONDA_ACTIVATION_CMD.replace('&','') 
        ALPHAFOLD_HOME = Plugin.getVar('ALPHAFOLD_HOME')
        ALPHAFOLD_DATABASE_DIR = Plugin.getVar('ALPHAFOLD_DATABASE_DIR')
        OUTPUT_DIR = os.path.abspath(self._getExtraPath())
        inputFastaFile = self.createInputFastaFile(seqs)
        if len(seqs) > 1:
            multimer = 'multimer'
        else:
            multimer = 'monomer'
        if isProkaryote:
            isProkaryote = '-l true'
        else:
            isProkaryote = ''
        if doGpu:
            gpu = "-a %s"  % gpusToUse
        else:
            gpu = ''
        command=f"""#!/bin/bash
{CONDA_ACTIVATION_CMD} # activate conda comamnd
conda activate alphafold # activate conda alphafold enviroment
# The following environment variable settings may help
# with larger polypeptide calculations (> 1,200 aa).
TF_FORCE_UNIFIED_MEMORY=1
XLA_PYTHON_CLIENT_MEM_FRACTION=0.5
XLA_PYTHON_CLIENT_ALLOCATOR=platform
#
cd {ALPHAFOLD_HOME}
/bin/bash {ALPHAFOLD_HOME}/run_alphafold.sh  \
-d {ALPHAFOLD_DATABASE_DIR} \
-o {OUTPUT_DIR} \
-f {inputFastaFile} \
-m {multimer} \
{isProkaryote} \
-t {maxTemplateDate} \
{gpu} \
{extraFlags}
"""
        alphaFoldScriptName = os.path.abspath(self._getExtraPath("alphafold.sh"))
        f = open(alphaFoldScriptName, "w")
        f.write(command)
        f.close()
        self.runJob('/bin/bash', alphaFoldScriptName, cwd=ALPHAFOLD_HOME) 

        outFileNames = []
        outputdir = self._getExtraPath(self.INPUTFASTAFILE)
        searchPattern = os.path.join(outputdir, "ranked_?.pdb")
        for outFileName in sorted(glob.glob(searchPattern)):
            outFileNames.append(outFileName)
        # check if we have any output
        if not outFileNames:
            error_message = "No atomic model selected"
            raise Exception(error_message)
        else:
            self.createOutputStep(outFileNames)


    def _getModelFromEBI(self, uniProtID):
        '''Fetch structures from EBI AlphaFold database using UniProt sequence ID.
           Example for UniProt P29474.
           https://alphafold.ebi.ac.uk/files/AF-P29474-F1-model_v1.cif'''
        
        # get alphafold EBI database url
        self._get_alphafold_database_settings()
        data = {'uniprot_id': uniProtID, 'version': self.settings['database_version']}
        # get model
        extraDir = self._getExtraPath()
        model_url = self.settings['database_url'].format(**data)
        outFileName = os.path.join(extraDir, uniProtID + ".cif")
        status, msg = fetch_file(model_url, retry = 3, json=False, outFilename=outFileName)
        if not status:
            error_message = f"ERROR: Can not retieve {uniProtID} from EBI Alphafold database. Is {uniProtID} a valid UNIPROT ID?"
            error_message += f"SYSTEM report error {msg}"
            raise Exception(error_message)
        else:
            self.createOutputStep([outFileName])

    def _get_alphafold_database_settings(self):
        """ get alphafold database settings from 
             https://www.rbvi.ucsf.edu/chimerax/data/status/alphafold_database.json"""
        url = 'https://www.rbvi.ucsf.edu/chimerax/data/status/alphafold_database.json'
        self.settings = {}
        try:
            self.settings = fetch_file(url, json=True, outFilename=None)
        except Exception:
            print("Could not reach update site")

        if not self.settings:
            raise Exception('No alphafold database settings found')

    def _getModelFromBlast(self, sequence_data, hideMessage):
        """run a blast and search model for 5 closest matches
        :param text sequence_data: sequence in fasta format

        We will use chimera for this since it is faster than NCBI and
        allows you to search in the alfafold database 
        """
        # create script chimera
        dim = 150  # eventually we will create a PDB library that
                   # computes PDB dim
        sampling = 1.
        tmpFileName = os.path.abspath(self._getTmpPath("axis_input.bild"))
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=tmpFileName,
                                         sampling=sampling)
        chimeraScriptFileName = "chimeraPythonScript.py"
        f = open(self._getTmpPath(chimeraScriptFileName), "w")
        f.write('from chimerax.core.commands import run\n')

        f.write("run(session, 'open %s')\n" % tmpFileName)
        f.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates
        f.write("run(session, 'alphafold search %s')\n" % sequence_data)
        # Show help window
        if not hideMessage:
            msg = """Select desired homologous sequence and save the
corresponding atomic model with the command
scipionwrite #modelID [prefix myprefix]"""
            f.write(f"""
session.logger.error('''{msg}''')
""")
        
        # run the script:
        _chimeraScriptFileName = os.path.abspath(
            self._getTmpPath(chimeraScriptFileName))
        if len(self.extraCommands.get()) > 2:
            # TODO: parse extra command
            f.write(self.extraCommands.get())
            args = " --nogui " + _chimeraScriptFileName
        else:
            args = " " + _chimeraScriptFileName
        f.close()

        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)

        # run in the background
        cwd = os.path.abspath(self._getExtraPath())
        Plugin.runChimeraProgram(Plugin.getProgram(), args, 
                                 cwd=cwd, extraEnv=getEnvDictionary(self))
        outFileNames = []
        searchPattern = self._getExtraPath("*.cif")
        for outFileName in glob.glob(searchPattern):
            outFileNames.append(outFileName)
        # check if we have any output 
        if not outFileNames:
            error_message = f"No atomic model selected"
            raise Exception(error_message)
        else:
            self.createOutputStep(outFileNames)

    def _getModelFromColab(self, sequence_data, colabID, hideMessage, useTemplatesFromPDB=-1):
        """run colab to get an alphafold prediction
        We will use chimera for this.
        """
        # connect to localfold, we need to create a QT browser
        # QT is available in chimera's python
        colabScriptFileName = os.path.abspath(self._getExtraPath("colab.py"))
        f = open(self._getTmpPath(colabScriptFileName), "w")
        injectJavaScriptList = []

        ###
        # 1 CASE
        # monomer, chimera
        ###
        if colabID == self.CHIMERA:
            injectJavaScriptList.append(
                f'''document.querySelector("paper-input").setAttribute("value", "{sequence_data}");  + 
                   document.querySelector("paper-input").dispatchEvent(new Event("change"));
                ''')
            injectJavaScriptList.append('document.querySelector("colab-run-button").click()')
        ###
        # 2 CASE 
        # phenix multimer, reuse result
        ###
        elif colabID == self.PHENIX:  
            #counter = 0
            objId = self.getObjId()
            injectJavaScriptList.append(
                f'''document.querySelector("paper-input.flex[aria-labelledby='formwidget-1-label']").setAttribute("value", "{sequence_data}");'''
            )
            injectJavaScriptList.append(
                f'''document.querySelector("paper-input").dispatchEvent(new Event("change"));'''
            )

            injectJavaScriptList.append(
                f'''document.querySelector("paper-input.flex[aria-labelledby='formwidget-2-label']").setAttribute("value", "{objId}");'''
            )
            injectJavaScriptList.append(            
                f'''document.querySelector("paper-input").dispatchEvent(new Event("change"));'''
            )
            #injectJavaScriptList.append(f'document.querySelectorAll("colab-run-button")[{counter}].click()')
            #counter += 1
            if useTemplatesFromPDB>0:
                injectJavaScriptList.append(            
                    '''document.querySelector("input[aria-labelledby=formwidget-5-label]").click()'''
                )
                injectJavaScriptList.append(            
                    f'''document.querySelector("paper-input.flex[aria-labelledby='formwidget-6-label']").setAttribute("value", "{useTemplatesFromPDB}")'''
                )
                injectJavaScriptList.append(            
                    '''document.querySelector("input[aria-labelledby=formwidget-7-label]").click()'''
                )

        elif colabID == self.TEST:
            injectJavaScriptList.append('document.querySelector("colab-run-button").click()')


        createcolabscript = createColabScript(scriptFilePointer=f,
                                              url=self.url[colabID],
                                              extraPath=os.path.abspath(self._getExtraPath()),
                                              injectJavaScriptList=injectJavaScriptList,
                                              )
        f.close()
        args = " --nogui " + colabScriptFileName
        cwd = os.path.abspath(self._getExtraPath())
        Plugin.runChimeraProgram(Plugin.getProgram(), args, 
                                 cwd=cwd, extraEnv=getEnvDictionary(self))


        # uncompress file should be in results.zip in the extra directory
        # create results dir
        # uncompress data in
        # best model should go to output
        # identify other different from best
        # other accesible though  chimera

        ########################################
        # create script chimera (move to visualize)
        # dim = 150  # eventually we will create a PDB library that
        #            # computes PDB dim
        # sampling = 1.
        # tmpFileName = os.path.abspath(self._getTmpPath("axis_input.bild"))
        # # chimeraScriptFileName = "chimeraScript.py"
        # Chimera.createCoordinateAxisFile(dim,
        #                                  bildFileName=tmpFileName,
        #                                  sampling=sampling)
        # chimeraScriptFileName = "chimeraPythonScript.py"
        # f = open(self._getTmpPath(chimeraScriptFileName), "w")
        # f.write('from chimerax.core.commands import run\n')

        # f.write("run(session, 'open %s')\n" % 'best_model.pdb')
        # f.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates
        # f.write("color bfactor palette alphafold\n")
        # f.write("key red:low orange: yellow: cornflowerblue: blue:high\n")
        #######################################################
#         if not hideMessage:
#             title = "help"
#             msg = """Wait untill google colab ends and then
# close chimera,
# Some complementary information is
#  avaialble at ~/Downloads/ChimeraX/Alphafold"""
#             f.write(f"""
# session.logger.error('''{msg}''')
# """)
        f.close()

#         self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)
        outFileNames = []
        bestModelFileName = self._getExtraPath(os.path.join('results', 'best_model.pdb'))
        outFileNames.append(bestModelFileName)

        if not outFileNames:
            error_message = f"No atomic model selected"
            raise Exception(error_message)
        else:
            self.createOutputStep(outFileNames)

    def createOutputStep(self, atomStructPaths):
        """ Copy the atomic  structure and register the output object.
        :param list_string atomStructPath: list of atom struct files to be
                                        saved
        """
        for atomStructPath in atomStructPaths:
            if not exists(atomStructPath):
                raise Exception("Atomic structure not found at *%s*" % atomStructPath)
            if atomStructPath.endswith(".pdb") or atomStructPath.endswith(".cif"):
                pdb = emobj.AtomStruct()
                pdb.setFileName(atomStructPath)
                atomStructPath = os.path.basename(atomStructPath)
                if atomStructPath.endswith(".cif"):
                    keyword = atomStructPath.split(".cif")[0].replace(".","_")
                else:
                    keyword = atomStructPath.split(".pdb")[0].replace(".", "_")
                kwargs = {keyword: pdb}
                self._defineOutputs(**kwargs)


    # --------------------------- INFO functions ---------------------------

    def _validate(self):
        errors = []
        if self.source == self.IMPORT_LOCAL_ALPHAFOLD:        
            ALPHAFOLD_HOME = Plugin.getVar('ALPHAFOLD_HOME')
            binary = os.path.join(ALPHAFOLD_HOME, "run_alphafold.sh") 
            if not os.path.exists(binary):
                errors.append("No valid Alphafold instalation found")
        return errors

    def _citations(self):
        return ['Alphafold2021']


    def _summary(self):
        summary = []
        summary.append('PLDDTs (monomer) or iptm+ptm (multimer)')
        try:
            outputdir = self._getExtraPath(self.INPUTFASTAFILE)
            fileName = os.path.join(outputdir, 'ranking_debug.json')
            if os.path.exists(fileName):
                self._log.info('inside if')
                f = open(fileName)
                data = json.load(f)
                for i, model in enumerate(data['order']):
                    # if multimer:
                        # label = 'iptm+ptm'
                    #else:
                         # label = 'pldots'
                    if 'plddts' in data:
                        key = 'plddts'
                    elif 'iptm+ptm' in data:
                        key = 'ptm+ptm'
                    summary.append("*%d* %.2f" % (i, data[key][model]))
            else:
                summary.append('alphafold ranking not yet computed')
        except:
            summary.append('Cannot create summary')
        return summary 

#----------------- utils -----------------


def _parse_uniprot_id(unitProtID):
    """Check uniprotid. Call from validate """
    if len(unitProtID) not in (6, 10):
        raise Exception('UniProt identifiers must be 6 or 10 characters long, got "%s"'
            % unitProtID)
    return unitProtID.upper()

def fetch_file(url, retry = 3, json=False, outFilename=None):
    """ fetch file from url, retry 'retries' number of times
    :param str url: full url to file
    :param int reties: number of attemps
    :param bool json: parse json response
    :param outFilename
    """
    for r in range(retry):
        try:
            n = os.path.basename(url)
            response = requests.get(url) # Downloading the file and saving it at app/test with the file name n
        except Exception as e:
            if r < 2:
                print(f'Failed. Attempt # {r + 1}')
            else:
                print('Error encoutered at third attempt downloading {n}')
                print(e)
        else:
            print(f"Success: {n}")
            break

    if json:
        return response.json()

    if outFilename is not None:
        # check if we have a valid answer
        if response.text.find('<?xml version=') != -1:
            # error case, parse return xml string
            root = ET.fromstring(response.text)
            errorMessage = root.find("Message")
            return False, errorMessage.text
        with open(outFilename, 'wb') as f:
            f.write(response.content)
            return True, ''
    else:  # return content
        return True, response.content

# get directory in which alphafold model will be downloaded
# this is weak since if two alphafold are run in parallel
# it will fail
def _unique_download_directory():
    from os.path import expanduser, join, exists
    ddir = expanduser('~/Downloads')
    adir = join(ddir, 'ChimeraX', 'AlphaFold')
    from os import makedirs
    makedirs(adir, exist_ok = True)
    for i in range(1,1000000):
        path = join(adir, 'prediction_%d' % i)
        if not exists(path):
            path = join(adir, 'prediction_%d' % (i-1))
            break
    return path

def _getSize(file_path):
    if os.path.isfile(file_path): 
        st = os.stat(file_path)
        return st.st_size
    else:
        return -1

def _waitForFile(file_path):
    counter = 0
    current_size = _getSize(file_path)
    time.sleep(60)
    while current_size !=_getSize(file_path) or _getSize(file_path)==0 or not os.path.exists(file_path):
        current_size =_getSize(file_path)
        time.sleep(60)# wait download
        counter += 1
        if counter > 240:  # break after four hours
            break

def _findDownloadDirAndGetModels():
    "Return last subdirectory created by alphafold-colab"
    import glob
    downloadDir = _unique_download_directory()
    return [os.path.join(downloadDir, 'best_model.pdb')] + glob.glob(os.path.join(downloadDir, 'model_?_unrelaxed.pdb'))
