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

from pyworkflow import VERSION_1_2
try:
    from pyworkflow.em.data import AtomStruct
except:
    from pyworkflow.em.data import PdbFile as AtomStruct

from pyworkflow.em import Volume
from pyworkflow.em.convert import ImageHandler
from pyworkflow.em.data import Transform
from pyworkflow.em.convert.headers import Ccp4Header
from pyworkflow.em.protocol import EMProtocol

from pyworkflow.em.viewers.viewer_chimera import (Chimera,
                                                  chimeraScriptFileName,
                                                  sessionFile,
                                                  chimeraMapTemplateFileName,
                                                  chimeraScriptFileName,
                                                  chimeraPdbTemplateFileName)

from pyworkflow.protocol.params import (MultiPointerParam,
                                        PointerParam,
                                        StringParam)
from pyworkflow.utils.properties import Message

from chimera import Plugin


class ChimeraProtBase(EMProtocol):
    """Protocol to perform rigid fit using Chimera.
        Execute command *scipionwrite [model #n] [refmodel #p]
        [saverefmodel 0|1]* from command line in order to transferm fitted
        pdb to scipion. Default values are model=#0,
        refmodel =#1 and saverefmodel 0 (false).
        model refers to the pdb file. refmodel to a 3Dmap"""
    _version = VERSION_1_2

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=True):
        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input Volume', allowsNull=True,
                      help="Volume to process")
        form.addParam('pdbFileToBeRefined', PointerParam,
                      pointerClass="AtomStruct",
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
            form.addLine('''Execute command *scipionwrite [model #n] [refmodel #p] 
            [saverefmodel 0|1]* from command line in order to transfer structures 
            and 3D map volumes to SCIPION. 
            In the particular case in which you have only a volume and a structure, 
            default values are model #2, refmodel #1 and saverefmodel 0 (false). 
            Model refers to the PDBx/mmCIF file, refmodel to a 3D map volume. 
            If you have several structures and no volumes, you can save 
            all of them by executing commands *scipionwrite [model #1]*, 
            *scipionwrite [model #2]*, *scipionwrite [model #3]*, and so on.
            When you use the command line scipionwrite, the Chimera session will 
            be saved by default. Additionally, you can save the Chimera session 
            whenever you want by executing the command *scipionss". You will be 
            able to restore the saved session by using the protocol chimera restore 
            session (SCIPION menu: Tools/Calculators/chimera restore session). ''')

        return form #DO NOT remove this return

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
                         self._getExtraPath(sessionFile)
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
            pdbModelCounter += 1
            x_input, y_input, z_input = _inputVol.getShiftsFromOrigin()
            inputVolFileName = os.path.abspath(ImageHandler.removeFileType(
                _inputVol.getFileName()))
            f.write("runCommand('open %s')\n" % inputVolFileName)
            f.write("runCommand('volume #1 style surface voxelSize %f')\n"
                    % _inputVol.getSamplingRate())
            f.write("runCommand('volume #1 origin %0.2f,%0.2f,%0.2f')\n"
                    % (x_input, y_input, z_input))

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

        # outvolName, this volume may or may not exists
        volFileName = self._getExtraPath((chimeraMapTemplateFileName) % 1)

        # if we have outvol
        if os.path.exists(volFileName):
            # if we do not have an explicit inputvol check if it
            # is a volume associated to the  pdb
            if self.inputVolume.get() is None:
                _inputVol = self.pdbFileToBeRefined.get().getVolume()
            else:
                _inputVol = self.inputVolume.get()
            if _inputVol is not None: # we have inputVol
                oldSampling = _inputVol.getSamplingRate()

            vol = Volume()  # this is an output volume object
            vol.setLocation(volFileName)

            # fix mrc header
            ccp4header = Ccp4Header(volFileName, readHeader=True)
            sampling = ccp4header.computeSampling()
            vol.setSamplingRate(sampling)

            #find origin
            if _inputVol is not None:
                if self.inputVolume.get() is None:
                    origin = self.pdbFileToBeRefined.get().getVolume(). \
                        getOrigin(force=True)
                else:
                    origin = self.inputVolume.get().getOrigin(
                        force=True)

                newOrigin = vol.originResampled(origin, oldSampling)
                vol.setOrigin(newOrigin)

            else: # in this case there is no inputvol
                  # but there is outputvol.
                if self.pdbFileToBeRefined.get().hasOrigin():
                    origin = self.pdbFileToBeRefined.get().getOrigin()
                else:
                    origin = Transform()
                    shifts = ccp4header.getOrigin()
                    origin.setShiftsTuple(shifts)

                if origin is None:
                    origin = vol.getOrigin(force=True)
                vol.setOrigin(origin)

            self._defineOutputs(output3Dmap=vol)

            if self.inputVolume.get() is None:
                self._defineSourceRelation(
                    self.pdbFileToBeRefined.get(), vol)
            else:
                self._defineSourceRelation(self.inputVolume.get(), vol)
        #we do not have output volume
        else:
            if self.inputVolume.get() is None:
                vol = self.pdbFileToBeRefined.get().getVolume()
            else:
                vol = self.inputVolume.get()

        # Now check pdb files
        directory = self._getExtraPath()
        counter = 1
        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".pdb") or filename.endswith(".cif"):
                path = os.path.join(directory, filename)
                pdb = AtomStruct()
                pdb.setFileName(path)
                if vol is not None:
                    pdb.setVolume(vol)
                keyword = "outputPdb_%02d" % counter
                counter += 1
                kwargs = {keyword: pdb}
                self._defineOutputs(**kwargs)
                self._defineSourceRelation(self.pdbFileToBeRefined.get(), pdb)
                for pdbFile in self.inputPdbFiles:
                    self._defineSourceRelation(pdbFile.get(), pdb)

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
        if (self.getOutputsSize() > 0):
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
        methodsMsgs = []
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
import chimera
from chimera import runCommand
def newFileName(template):
    counter = 1
    while os.path.isfile(template%counter):
        counter += 1
    return template%counter


def saveSession(sessionFileName):
    runCommand('save %s' % sessionFileName)

def restoreSession(sessionFileName):
    runCommand('open %s' % sessionFileName)

def saveModel(model, refModel, fileName):
    runCommand('write relative %s %s %s'%(refModel, model, fileName))
    
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
scipionwrite model #n refmodel #p saverefmodel 0/1
"""
def cmd_scipionWrite(scipionWrite,args):
  from Midas.midas_text import doExtensionFunc
'''

chimeraScriptMain = '''
  def scipionWrite(model="#%(pdbID)d",refmodel="#%(_3DmapId)d",
     saverefmodel=0):
     # get model (pdb) id
     try:
         saveSession('%(sessionFileName)s')
     except Exception as e:
         f = open ('/tmp/chimera_error.txt','w')
         f.write(e.message)
         f.close()
         
     # modelId = int(model[1:])# model to write  1
     refModelId = int(refmodel[1:])# coordenate system refers to this model 0

     # get actual models
     # model    = chimera.openModels.list()[modelId]
     # TODO: this ID is wrong if models before refmodel are modified
     refModel = chimera.openModels.list()[refModelId]

     #for m in chimera.openModels.list():
     #    if m.id != modelId:
     #        continue
     #    else:
     #        break
     # Save the PDB relative to the volume coordinate system
     # TODO: check if this Will work if the reference is a PDB?
     # from Midas import write
     fileName = newFileName('%(pdbFileTemplate)s')
     saveModel(model, refModel, fileName)
     # alternative way to save  the pdb file using a command
     # run('write relative #1 #0 pdb_path')

     # Save the map if sampling rate has been changed
     if saverefmodel:
         _save_grid_data(refModel.data,"%(chimeraMapTemplateFileName)s" )
     # always save session when write
     # beep(0.1)

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
addCommand('scipionwrite', cmd_scipionWrite, help="http://scipion.cnb.csic.es")
addCommand('scipionss', cmd_scipionSaveSession, help="http://scipion.cnb.csic.es")
addCommand('scipionrs', cmd_scipionRestoreSession, help="http://scipion.cnb.csic.es")
'''


def createScriptFile(pdbID, _3DmapId,
                     pdbFileTemplate, mapFileTemplate,
                     f, sessionFileName=''):
    f.write(chimeraScriptHeader)
    d = {}
    d['pdbID'] = pdbID
    d['_3DmapId'] = _3DmapId
    d['pdbFileTemplate'] = pdbFileTemplate  # % 1
    d['chimeraMapTemplateFileName'] = mapFileTemplate % 1
    d['sessionFileName'] = sessionFileName
    f.write(chimeraScriptMain % d)
