# **************************************************************************
# *
# * Authors:     Roberto Marabini (roberto@cnb.csic.es)
# *              Yunior C. Fonseca Reyna (cfonseca@cnb.csic.es)
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

import pwem
import pyworkflow.utils as pwutils

from .constants import CHIMERA_HOME, CHIMERA_HEADLESS_HOME, V1_13_1

__version__ = "3.0.1"
_logo = "chimera_logo.png"
_references = ['Pettersen2004']


class Plugin(pwem.Plugin):
    _homeVar = CHIMERA_HOME
    _pathVars = [CHIMERA_HOME]
    _supportedVersions = V1_13_1

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(CHIMERA_HOME, 'chimera-1.13.1')
        cls._defineEmVar(CHIMERA_HEADLESS_HOME, 'chimera_headless')

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        d = {}
        d['PATH'] = cls.getHome('bin')
        if "REMOTE_MESA_LIB" in os.environ:
            d["LD_LIBRARY_PATH"] = os.environ['REMOTE_MESA_LIB']
        environ.update(d, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def runChimeraProgram(cls, program, args="", cwd=None):
        """ Internal shortcut function to launch chimera program. """
        env = cls.getEnviron()
        pwutils.runJob(None, program, args, env=env, cwd=cwd)

    @classmethod
    def getProgram(cls, progName="chimera"):
        """ Return the program binary that will be used. """
        cmd = cls.getHome('bin', progName)
        return str(cmd)

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(V1_13_1)

    @classmethod
    def defineBinaries(cls, env):

        SW_CH = env.getEmFolder()
        chimera_1_13_1_command = [('./scipion_installer',
                            '%s/chimera-1.13.1/bin/chimera' % SW_CH)]

        env.addPackage('chimera', version='1.13.1',
                       tar='chimera-1.13.1-linux_x86_64.tgz',
                       commands=chimera_1_13_1_command,
                       default=True)
