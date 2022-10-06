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
from ..utils import getEnvDictionary

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
import configparser
import shutil

class ChimeraProtBase(EMProtocol):
    """Base class  for chimera protocol"""
    _version = VERSION_3_0
    @classmethod
    def getClassPackageName(cls):
        return "chimerax"

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
            # form.addLine(''' scipionwrite model #n [refmodel #p] [prefix stringAddedToFilename]
            form.addLine(''' scipionwrite #n [prefix stringAddedToFilename]
            scipionss
            scipionrs
            scipioncombine #n1,n2,n3... [modelid StringArg]
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
        _inputVol = None
        if (hasattr(self, 'inputVolume') and
                (self.inputVolume.get() is None) and
                (self.pdbFileToBeRefined.get() is not None)):
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
        f.write("open %s\n" % tmpFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates

        # input vol with its origin coordinates
        pdbModelCounter = 1
        if _inputVol is not None:
            pdbModelCounter += 1
            x_input, y_input, z_input = _inputVol.getShiftsFromOrigin()
            inputVolFileName = os.path.abspath(ImageHandler.removeFileType(
                _inputVol.getFileName()))
            f.write("open %s\n" % inputVolFileName)
            f.write("volume #%d style surface voxelSize %f\n"
                    % (pdbModelCounter, _inputVol.getSamplingRate()))
            f.write("volume #%d origin %0.2f,%0.2f,%0.2f\n"
                    % (pdbModelCounter, x_input, y_input, z_input))

        if hasattr(self, 'inputVolumes') and \
                self.inputVolumes is not None:
            for vol in self.inputVolumes:
                pdbModelCounter += 1
                f.write("open %s\n" % os.path.abspath(vol.get().getFileName()))
                x, y, z = vol.get().getShiftsFromOrigin()
                f.write("volume #%d style surface voxelSize %f\n"
                        % (pdbModelCounter, vol.get().getSamplingRate()))
                f.write("volume #%d origin %0.2f,%0.2f,%0.2f\n"
                        % (pdbModelCounter, x, y, z))

        if self.pdbFileToBeRefined.get() is not None:
            pdbModelCounter += 1
            pdbFileToBeRefined = self.pdbFileToBeRefined.get()
            f.write("open %s\n" % os.path.abspath(
                pdbFileToBeRefined.getFileName()))
            if pdbFileToBeRefined.hasOrigin():
                x, y, z = (pdbFileToBeRefined.getOrigin().getShifts())
                f.write("move %0.2f,%0.2f,%0.2f model #%d "
                        "coord #0\n" % (x, y, z, pdbModelCounter))

        # other pdb files
        if hasattr(self, 'inputPdbFiles'):
            for pdb in self.inputPdbFiles:
                pdbModelCounter += 1
                f.write("open %s\n" % os.path.abspath(pdb.get(
                ).getFileName()))
                if pdb.get().hasOrigin():
                    x, y, z = pdb.get().getOrigin().getShifts()
                    f.write("move %0.2f,%0.2f,%0.2f model #%d "
                            "coord #0\n" % (x, y, z, pdbModelCounter))

        # run the text:
        _chimeraScriptFileName = os.path.abspath(
            self._getTmpPath(chimeraScriptFileName))
        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
            args = " --nogui " + _chimeraScriptFileName
        else:
            args = " " + _chimeraScriptFileName

        f.close()

        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)

        # run in the background
        cwd = os.path.abspath(self._getExtraPath())
        Plugin.runChimeraProgram(Plugin.getProgram(), args, cwd=cwd, extraEnv=getEnvDictionary(self))

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        # Check vol and pdb files
        directory = self._getExtraPath()
        for filename in sorted(os.listdir(directory)):
            if not filename.startswith("tmp"):
                # files starting with "tmp" will not be converted in scipion objects
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
                    if filename.endswith(".cif"):
                        keyword = filename.split(".cif")[0].replace(".","_")
                    else:
                        keyword = filename.split(".pdb")[0].replace(".", "_")
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
                errors.append("CHIMERA_HOME = %s" % Plugin.getHome())

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
        return ['Goddard2018']

    def is_tool(self, name):
        """Check whether `name` is on PATH."""
        from shutil import which
        return which(name) is not None
