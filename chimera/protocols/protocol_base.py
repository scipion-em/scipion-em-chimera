# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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

import os

from pyworkflow import VERSION_3_0

try:
    from pwem.objects import AtomStruct
except ImportError:
    from pwem.objects import PdbFile as AtomStruct

from pwem.objects import Volume
from pwem.emlib.image import ImageHandler
from pwem.objects import Transform
from pwem.convert.headers import Ccp4Header
from pwem.protocols import EMProtocol

from pwem.viewers.viewer_chimera import (Chimera,
                                         sessionFile,
                                         chimeraMapTemplateFileName,
                                         chimeraScriptFileName,
                                         chimeraPdbTemplateFileName)

from pyworkflow.protocol.params import (MultiPointerParam,
                                        PointerParam,
                                        StringParam)
from pyworkflow.utils.properties import Message

from .. import Plugin


class ChimeraProtBase(EMProtocol):
    """Base class  for chimera protocol"""
    _version = VERSION_3_0

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=True):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      important=True,
                      help="Volume to process")
        form.addParam('inputVolumes', MultiPointerParam, pointerClass="Volume",
                      label='Input additional Volumes', allowsNull=True,
                      help="Other Volumes")
        form.addParam('pdbFileToBeRefined', PointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      important=True,
                      label='Atomic structure',
                      help="PDBx/mmCIF file that you can save after operating "
                           "with it.")
        form.addParam('inputPdbFiles', MultiPointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      label='Other atomic structures',
                      help="In case you need to load more PDBx/mmCIF files, "
                           "you can load them here and save them after "
                           "operating with them.")
        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="Add extra commands in cmd file. Use for testing")
        if doHelp:
            form.addSection(label='Help')
            form.addLine(''' scipionwrite model #n [refmodel #p] [prefix stringAddedToFilename]
            scipionss
            scipionrs
            Type 'help command' in chimera command line for details (command is the command name)''')

        return form  # DO NOT remove this return

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):

        self._insertFunctionStep('prerequisitesStep')
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions ---------------------------

    def prerequisitesStep(self):
        """
        """
        pass

    def runChimeraStep(self):
        # building script file including the coordinate axes and the input
        # volume with samplingRate and Origin information
        f = open(self._getTmpPath(chimeraScriptFileName), "w")
        f.write("from chimera import runCommand\n")

        # create coherent header
        createScriptFile(1,  # model id pdb
                         1,  # model id 3D map
                         self._getExtraPath(chimeraPdbTemplateFileName),
                         self._getExtraPath(chimeraMapTemplateFileName),
                         f,
                         self._getExtraPath(sessionFile),
                         )

        if self.inputVolume.get() is None:
            _inputVol = self.pdbFileToBeRefined.get().getVolume()
        else:
            _inputVol = self.inputVolume.get()

        # building coordinate axes
        if _inputVol is not None:
            dim = _inputVol.getDim()[0]
            sampling = _inputVol.getSamplingRate()
        else:
            dim = 150  # eventually we will create a PDB library that
            # computes PDB dim
            sampling = 1.

        tmpFileName = os.path.abspath(self._getTmpPath("axis_input.bild"))
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=tmpFileName,
                                         sampling=sampling)
        f.write("runCommand('open %s')\n" % tmpFileName)
        f.write("runCommand('cofr 0,0,0')\n")  # set center of coordinates

        # input vol with its origin coordinates
        pdbModelCounter = 1
        if _inputVol is not None:
            x_input, y_input, z_input = _inputVol.getShiftsFromOrigin()
            inputVolFileName = os.path.abspath(ImageHandler.removeFileType(
                _inputVol.getFileName()))
            f.write("runCommand('open %s')\n" % inputVolFileName)
            f.write("runCommand('volume #%d style surface voxelSize %f')\n"
                    % (pdbModelCounter, _inputVol.getSamplingRate()))
            f.write("runCommand('volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                    % (pdbModelCounter, x_input, y_input, z_input))

        if self.inputVolumes is not None:
            if _inputVol is not None:
                pdbModelCounter += 1
            for vol in self.inputVolumes:
                f.write("runCommand('open %s')\n" % vol.get().getFileName())
                x, y, z = vol.get().getShiftsFromOrigin()
                f.write("runCommand('volume #%d style surface voxelSize %f')\n"
                        % (pdbModelCounter, vol.get().getSamplingRate()))
                f.write("runCommand('volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                        % (pdbModelCounter, x, y, z))
                pdbModelCounter += 1

        if self.pdbFileToBeRefined.get() is not None:
            pdbFileToBeRefined = self.pdbFileToBeRefined.get()
            f.write("runCommand('open %s')\n" % os.path.abspath(
                pdbFileToBeRefined.getFileName()))
            if pdbFileToBeRefined.hasOrigin():
                x, y, z = (pdbFileToBeRefined.getOrigin().getShifts())
                f.write("runCommand('move %0.2f,%0.2f,%0.2f model #%d "
                        "coord #0')\n" % (x, y, z, pdbModelCounter))

        # Alignment of sequence and structure
        if (hasattr(self, 'inputSequence') and
                hasattr(self, 'inputStructureChain')):
            if (self.inputSequence.get() is not None and
                    self.inputStructureChain.get() is not None):
                models = self.structureHandler.getModelsChains()
                if len(models) > 1:
                    f.write("runCommand('select #%d.%s:.%s')\n"
                            % (pdbModelCounter, str(self.selectedModel),
                               str(self.selectedChain)))
                else:
                    f.write("runCommand('select #%d:.%s')\n"
                            % (pdbModelCounter, str(self.selectedChain)))
                # f.write("runCommand('sequence selection')\n") # To open
                # only the
                # sequence of
                # the selected
                # chain
                if self._getOutFastaSequencesFile is not None:
                    alignmentFile = self._getOutFastaSequencesFile()
                    f.write("runCommand('open %s')\n" % alignmentFile)

        # other pdb files
        pdbModelCounter += 1
        for pdb in self.inputPdbFiles:
            f.write("runCommand('open %s')\n" % os.path.abspath(pdb.get(
            ).getFileName()))
            if pdb.get().hasOrigin():
                x, y, z = pdb.get().getOrigin().getShifts()
                f.write("runCommand('move %0.2f,%0.2f,%0.2f model #%d "
                        "coord #0')\n" % (x, y, z, pdbModelCounter))
            # TODO: Check this this this this this this
            pdbModelCounter += 1

        # run the text:
        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
            args = " --nogui --script " + self._getTmpPath(
                chimeraScriptFileName)
        else:
            args = " --script " + self._getTmpPath(chimeraScriptFileName)

        f.close()

        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)

        # run in the background
        Chimera.runProgram(Plugin.getProgram(), args)

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        # Check vol and pdb files
        directory = self._getExtraPath()
        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".mrc"):
                volFileName = os.path.join(directory, filename)
                vol = Volume()
                vol.setFileName(volFileName)

                # fix mrc header
                ccp4header = Ccp4Header(volFileName, readHeader=True)
                sampling = ccp4header.computeSampling()
                origin = Transform()
                shifts = ccp4header.getOrigin()
                origin.setShiftsTuple(shifts)
                vol.setOrigin(origin)
                vol.setSamplingRate(sampling)
                keyword = filename.split(".mrc")[0]
                kwargs = {keyword: vol}
                self._defineOutputs(**kwargs)

            if filename.endswith(".pdb") or filename.endswith(".cif"):
                path = os.path.join(directory, filename)
                pdb = AtomStruct()
                pdb.setFileName(path)
                keyword = filename.split(".pdb")[0].replace(".","_")
                kwargs = {keyword: pdb}
                self._defineOutputs(**kwargs)


    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = Plugin.getProgram()
        if program is None:
            errors.append("Missing variable CHIMERA_HOME")
        elif not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)

        # If there is any error at this point it is related to config variables
        if errors:
            errors.append("Check configuration file: ~/.config/scipion/"
                          "scipion.conf")
            errors.append("and set CHIMERA_HOME variables properly.")
            if program is not None:
                errors.append("Current value:")
                errors.append("CHIMERA_HOME = %s" % os.environ['CHIMERA_HOME'])

        return errors

    def _summary(self):
        # Think on how to update this summary with created PDB
        summary = []
        if self.getOutputsSize() > 0:
            directory = self._getExtraPath()
            counter = 1
            summary.append("Produced files:")
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".pdb"):
                    summary.append(filename)
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".mrc"):
                    summary.append(filename)
            summary.append("we have some result")
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        return summary

    def _methods(self):
        methodsMsgs = list()
        methodsMsgs.append("TODO")

        return methodsMsgs

    def _citations(self):
        return ['Pettersen2004']

    def is_tool(self, name):
        """Check whether `name` is on PATH."""
        from distutils.spawn import find_executable
        return find_executable(name) is not None


# define scipion_write command
chimeraScriptHeader = '''
import os
import ntpath
#import chimera
from chimera import runCommand as rc
from chimera import openModels
from chimera import NonChimeraError

def newFileName(template, model):
    # counter = 1
    # while os.path.isfile(template%counter):
    #    counter += 1
    #f = open('/tmp/kk', 'w')
    #f.write(template)
    #f.close()
    
    return template % model


def saveSession(sessionFileName):
    rc('save %s' % sessionFileName)

def restoreSession(sessionFileName):
    rc('open %s' % sessionFileName)

def saveModel(_model, _refModel, fileName, mapFileName, prefix):
    import _molecule
    import VolumeViewer
    if prefix is not None:
        prefix = prefix.replace(".", "_dot_")
    # parse modelID #id.subid
    modelId = int(_model.split('.')[0][1:]) # object to be saved. Either pdb or 3d map.
    modelSubId = 0
    if _model.find('.') != -1:
        modelSubId = int(_model[1]) 

    # get model so we can asses if it is a 3D map or an atomic structure
    model = openModels.list(id=modelId, subid=modelSubId)[0]
    
    if isinstance(model, _molecule.Molecule):
        if prefix is not None:
            fileName = os.path.join(ntpath.dirname(fileName),
                                    prefix + ntpath.basename(fileName))
        rc('write relative %s %s %s'%(_refModel, _model, fileName))
    elif isinstance(model, VolumeViewer.volume.Volume):
        if prefix is not None:
            mapFileName = os.path.join(ntpath.dirname(mapFileName),
                                       prefix + ntpath.basename(mapFileName))
        _save_grid_data(model.data,"%s" % mapFileName )
    else:
        raise NonChimeraError("Can not write model %s" % (_model))

    
def _save_grid_data(refModel_data, fileName):
    from VolumeData import save_grid_data
    save_grid_data(refModel_data, os.path.abspath(fileName))
    
def beep(time):
   """I simply do not know how to create a portable beep sound.
      This system call seems to work pretty well if you have sox
      installed"""
   try:
      command = "play --no-show-progress -n synth %f sin 880"%time
      os.system(command)
   except:
      pass

"""load the script in chimera for example chimera --script ScipioChimeraExt/
ChimeraExtension.py
from chimera command line type one of the following three options
scipionwrite
scipionwrite model #n
scipionwrite model #n refmodel #p
"""
def cmd_scipionWrite(scipionWrite,args):
  from Midas.midas_text import doExtensionFunc
'''

chimeraScriptMain = '''
  def scipionWrite(model="%(modelID)s",refModel="%(refModelID)s", prefix=None, savesession='1'):
     #  savesession 1 -> save session
     #  savesession 0 -> do not save session    
     # Save the PDB relative to the volume coordinate system
     # TODO: check if this Will work if the reference is a PDB?
     atomStructFileName = newFileName('%(pdbFileTemplate)s', model[1:])
     mapFileName = newFileName('%(chimeraMapTemplateFileName)s', model[1:])
     saveModel(model, refModel, atomStructFileName, mapFileName, prefix)
     # alternative way to save  the pdb file using a command
     
     if savesession == '1':
        savesession = True
     else:
        savesession = False
     try:
         if savessesion:
             if prefix is None or prefix != 'DONOTSAVESESSION':  # this is needed for the tests
                 saveSession('%(sessionFileName)s')
     except Exception as e:
         f = open ('/tmp/chimera_error.txt','w')
         f.write(e.message)
         f.close()

  doExtensionFunc(scipionWrite,args)

def cmd_scipionSaveSession(scipionSaveSession,args):
  from Midas.midas_text import doExtensionFunc
  def scipionSaveSession():
     saveSession('%(sessionFileName)s')
     beep(0.1)

  doExtensionFunc(scipionSaveSession,args)

def cmd_scipionRestoreSession(scipionRestoreSession,args):
  from Midas.midas_text import doExtensionFunc
  def scipionRestoreSession():
     restoreSession('%(sessionFileName)s')
     beep(0.1)

  doExtensionFunc(scipionRestoreSession,args)

from Midas.midas_text import addCommand
addCommand('scipionwrite', cmd_scipionWrite, help="https://github.com/scipion-em/scipion-em-chimera/wiki/scipionwrite")
addCommand('scipionw', cmd_scipionWrite, help="https://github.com/scipion-em/scipion-em-chimera/wiki/scipionwrite")
addCommand('scipionss', cmd_scipionSaveSession, help="https://github.com/scipion-em/scipion-em-chimera/wiki/scipionss")
addCommand('scipionrs', cmd_scipionRestoreSession, help="https://github.com/scipion-em/scipion-em-chimera/wiki/scipionrs")
'''


def createScriptFile(_modelID, _refModelID,
                     pdbFileTemplate, mapFileTemplate,
                     f, sessionFileName=''):
    f.write(chimeraScriptHeader)
    d = dict()
    d['modelID'] = "#%d" % _modelID
    d['refModelID'] = "#%d" % _refModelID
    d['pdbFileTemplate'] = pdbFileTemplate  # % 1
    d['chimeraMapTemplateFileName'] = mapFileTemplate # % 1
    d['sessionFileName'] = sessionFileName
    f.write(chimeraScriptMain % d)
