# **************************************************************************
# *
# * Authors:     Grigory Sharov (sharov@igbmc.fr)
# *              Marta Martinez (mmmtnez@cnb.csic.es)
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

from .protocol_base import ChimeraProtBase
from pyworkflow.protocol.params import (PathParam)
from chimera.utils import getEnvDictionary
import os
from .. import Plugin


class ChimeraProtOperate(ChimeraProtBase):
    # """This protocol provides access to Chimera and allows to save
    #     the result in Scipion framework.
    #     Execute command *scipionwrite [model #n] [refmodel #p]
    #     [saverefmodel 0|1]* from command line in order to transferm fitted
    #     pdb to scipion. Default values are model=#0,
    #     refmodel =#1 and saverefmodel 0 (false).
    #     model refers to the pdb file. refmodel to a 3Dmap"""
    """This protocol provides access to Chimera and allows to save the result
            in Scipion framework.
            Execute command *scipionwrite #n [prefix stringAddedToFilename]
            model refers to the pdb file"""

    """
        Chimera Operate (ChimeraProtOperate) — User Manual

        Overview

        The Chimera Operate protocol provides a flexible interface between
        Scipion and ChimeraX, allowing users to execute custom Chimera scripts
        within a Scipion workflow. It is designed for advanced structural
        visualization, manipulation, and automation tasks that require more
        flexibility than predefined protocols.

        Inputs and Workflow

        The protocol accepts atomic structures and volume maps inherited from
        the Chimera base protocol, together with an optional ChimeraX or Python
        script. If no script is provided, the protocol executes the default
        Chimera workflow. Otherwise, it loads the user script, replaces dynamic
        placeholders with the corresponding input file paths, generates a
        temporary execution script, and launches ChimeraX automatically.

        Script Execution

        The placeholders '{inputVolume}' and '{pdbFileToBeRefined}' are
        automatically replaced by the input volume and atomic structure paths,
        allowing reusable and dataset-independent scripts.

        Python scripts ending in '.py' are executed using ChimeraX script mode,
        while other files are interpreted as standard Chimera command scripts.
        Python execution enables advanced automation, conditional logic, and
        large-scale structural processing workflows.

        Outputs and Applications

        The protocol does not define fixed outputs because results depend on
        the operations implemented in the provided script. Users may generate
        modified atomic models, transformed density maps, Chimera sessions,
        rendered images, or custom structural analyses.

        Final Perspective

        Chimera Operate acts as a programmable bridge between Scipion and
        ChimeraX, enabling customized and reproducible structural biology
        workflows. Its flexibility makes it especially useful for advanced
        cryo-EM visualization, automated scripting, and specialized molecular
        analysis tasks.
        """
    _label = 'operate'

    def _defineParams(self, form):
        super(ChimeraProtOperate, self)._defineParams(form)
        section = form.getSection('Input')

        section.addParam(
            'pythonScript', PathParam,
            default='',
            label='Execute script in chimera',
            help="""Add a Python script to be executed in chimera.
'{inputVolume}' will be replaced by the input volume file name
and '{pdbFileToBeRefined}' by the pdb file name.
If script ends in "py" it will be executed as "chimerax --script scriptname"
else as "chimerax scriptname".
""")

    def runChimeraStep(self):
        self._log.info('Running Chimera Operate Protocol')
        chimeraScriptFileName = self.pythonScript.get().strip()
        if len(chimeraScriptFileName) < 3:
            super().runChimeraStep()
        else:
            # get pdb file name and input volume file name
            if self.pdbFileToBeRefined.get() is not None:
                pdbFileToBeRefined = os.path.abspath(
                    self.pdbFileToBeRefined.get().getFileName())
            else:
                pdbFileToBeRefined = ''
            if self.inputVolume.get() is not None:
                inputVolume = os.path.abspath(
                    self.inputVolume.get().getFileName())
            else:
                inputVolume = ''
            # open script file,read content and replace placeholders
            with open(chimeraScriptFileName, "r") as f:
                script = f.read()
                script = script.replace(
                    '{inputVolume}', inputVolume)
                script = script.replace(
                    '{pdbFileToBeRefined}', pdbFileToBeRefined)
            # create a temporary script file with the modified content
            _chimeraScriptFileName = os.path.basename(
                chimeraScriptFileName)  # just the file name
            tmpChimeraScriptFileName = os.path.abspath(self._getExtraPath(
                _chimeraScriptFileName))
            with open(tmpChimeraScriptFileName, "w") as f:
                f.write(script)

            # run chimera with the temporary script file
            if tmpChimeraScriptFileName[-2:] == 'py':
                args = " --script " + tmpChimeraScriptFileName
            else:
                args = " " + tmpChimeraScriptFileName
            self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)

            # run in the background
            cwd = os.path.abspath(self._getExtraPath())
            Plugin.runChimeraProgram(
                Plugin.getProgram(), args, cwd=cwd,
                extraEnv=getEnvDictionary(self))
