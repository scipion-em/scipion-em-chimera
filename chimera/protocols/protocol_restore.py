# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pyworkflow.protocol.params import PointerParam, StringParam
import os
from pwem.viewers.viewer_chimera import sessionFile
from .protocol_base import ChimeraProtBase
from chimera import Plugin

from chimera.utils import getEnvDictionary


class ChimeraProtRestore(ChimeraProtBase):
    """This protocol opens Chimera and restores a session
      that has been stored each time a 3Dmap or an atomic structure 
      by using `scipionwrite` or `scipionss` commad.
      Execute command *scipionwrite #n [prefix stringAddedToFilename]
      model refers to the pdb file"""

    """
        Chimera Restore Session (ChimeraProtRestore) — User Manual

        Overview

        The Chimera Restore Session protocol allows users to reopen and
        continue previously saved ChimeraX sessions directly inside the
        Scipion framework. Its main purpose is to recover visualization,
        structural analysis, and modeling environments that were stored
        during earlier Chimera-based workflows using the `scipionwrite`
        or `scipionss` commands.

        In cryo-EM and structural biology workflows, this protocol is
        especially useful for continuing interactive analysis sessions,
        recovering molecular visualization states, or revisiting complex
        fitting and modeling procedures without rebuilding the workspace
        manually.

        Inputs and Workflow

        The protocol requires a previous Chimera-related protocol as input,
        including rigid fitting, operate, modeling, subtraction, or previous
        restore sessions. During execution, the protocol retrieves the saved
        Chimera session associated with the selected protocol and restores
        the corresponding molecular and volumetric environment.

        The workflow automatically inherits relevant inputs such as atomic
        structures, reference maps, and additional PDB files from the parent
        protocol. This ensures consistency between the original session and
        the restored environment.

        Session Restoration

        The protocol locates the Chimera session file generated during the
        original execution and launches ChimeraX using that stored session.
        As a result, visualization parameters, loaded models, map positions,
        orientations, and user modifications are recovered automatically.

        This functionality is particularly valuable in long or interactive
        structural biology workflows where manual reconstruction of the
        analysis environment would be time-consuming or error-prone.

        Validation and Consistency Checks

        Before execution, the protocol verifies that a valid session file
        exists in the selected parent protocol directory. If no session has
        been previously saved, execution is stopped and an error is reported.

        This validation step prevents invalid restoration attempts and helps
        maintain workflow reproducibility within Scipion projects.

        Final Perspective

        Chimera Restore Session acts as a session recovery and continuity
        mechanism between ChimeraX and Scipion workflows. Rather than
        performing new structural calculations, it provides a reproducible
        way to reopen and continue previous molecular visualization and
        analysis environments, making it especially useful in iterative
        cryo-EM interpretation and structural modeling workflows.
        """
        # Execute command *scipionwrite [model #n] [refmodel #p]
        # [saverefmodel 0|1]* from command line in order to transfer fitted
        # pdb to scipion. Default values are model=#0,
        # refmodel =#1 and saverefmodel 0 (false).
        # model refers to the pdb file. refmodel to a 3Dmap"""

    _label = 'restore session'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputProtocol', PointerParam,
                      label="Input protocols", important=True,
                      pointerClass='ChimeraProtOperate, '
                                   'ChimeraProtRigidFit, '
                                   'ChimeraModelFromTemplate, '
                                   'ChimeraProtRestore, '
                                   'ChimeraSubtractionMaps',
                      help="protocol to be reloaded")

        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="Add extra commands in cmd file. Use for testing")
        form.addSection(label='Help')
        form.addLine('''Execute command *scipionwrite [model #n] [refmodel #p]
        [saverefmodel 0|1]* from command line in order to transfer structures
        and 3D map volumes to SCIPION.
        In the particular case in which you have only a volume and a structure,
        default values are model #2, refmodel #1 and saverefmodel 0 (false).
        Model refers to the PDBx/mmCIF file, refmodel to a 3D map volume.
        If you have several structures and no volumes, you can save
        all of them by executing commands *scipionwrite [model #1]*,
        *scipionwrite [model #2]*, *scipionwrite [model #3]*, and so on.''')

    # --------------------------- INSERT steps functions --------------------

    def prerequisitesStep(self):
        """
        """
        self.parentProt = self.inputProtocol.get()
        self.parentProt.setProject(self.getProject())  # I do not really
        # understand this line

        if hasattr(self, 'inputVolume') and \
                self.inputVolume is not None:
            self.inputVolume = self.parentProt.inputVolume
        if hasattr(self, 'pdbFileToBeRefined') and \
                self.pdbFileToBeRefined is not None:
            self.pdbFileToBeRefined = self.parentProt.pdbFileToBeRefined
        if hasattr(self, 'inputPdbFiles') and \
                self.inputPdbFiles is not None:
            self.inputPdbFiles = self.parentProt.inputPdbFiles

    def runChimeraStep(self):
        # create CMD file
        parentSessionFileName = self.parentProt._getExtraPath(sessionFile)

        # if len(self.extraCommands.get()) > 2:
        #     f.write(self.extraCommands.get())
        #     args = " --nogui --cmd " + self._getTmpPath(
        #         chimeraScriptFileName)

        # run in the background
        cwd = os.path.abspath(self._getExtraPath())
        Plugin.runChimeraProgram(Plugin.getProgram(), os.path.abspath(parentSessionFileName), cwd=cwd, extraEnv=getEnvDictionary(self))

    def _validate(self):
        errors = super(ChimeraProtRestore, self)._validate()
        parentProt = self.inputProtocol.get()
        parentProt.setProject(self.getProject())  # I do not really understand
        # this line
        sessionFileName = parentProt._getExtraPath(sessionFile)
        # Check SESSION.py exists
        if not os.path.exists(sessionFileName):
            errors.append("Error: No session saved by protocol: %s\n"
                          % parentProt.getObjLabel())

        return errors
