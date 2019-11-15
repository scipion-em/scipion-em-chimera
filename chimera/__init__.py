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

import pyworkflow.em
import pyworkflow.utils as pwutils


from bibtex import _bibtex # Load bibtex dict with references

from chimera.constants import CHIMERA_HOME, CHIMERA_HEADLESS_HOME, V1_10_1, MAXIT_HOME, MAXIT_TAR, MAXIT_URL, MAXIT

_logo = "chimera_logo.png"

_references = ['Pettersen2004']

# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented


class Plugin(pyworkflow.em.Plugin):
    _homeVar = CHIMERA_HOME
    _pathVars = [CHIMERA_HOME]
    _supportedVersions = V1_10_1

    @classmethod
    def getMaxitHome(cls):
        return cls.getVar(MAXIT_HOME)

    @classmethod
    def getMaxitBin(cls):
        return os.path.join(cls.getMaxitHome(), 'bin', MAXIT)

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(CHIMERA_HOME, 'chimera-1.13.1')
        cls._defineEmVar(CHIMERA_HEADLESS_HOME, 'chimera_headless')
        cls._defineEmVar(MAXIT_HOME, os.path.join('maxit-10.1'))
                                                  # 'maxit-10.1'))

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        environ.update({'PATH': cls.getHome('bin'),
                        'LD_LIBRARY_PATH': os.environ['REMOTE_MESA_LIB'],
                        }, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def runChimeraProgram(cls, program, args=""):
        """ Internal shortcut function to launch chimera program. """
        env = cls.getEnviron()
        pwutils.runJob(None, program, args, env=env)

    @classmethod
    def getProgram(cls, progName="chimera"):
        """ Return the program binary that will be used. """
        cmd = cls.getHome('bin', progName)
        return str(cmd)

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(V1_10_1)

    @classmethod
    def defineBinaries(cls, env):

        SW_CH = env.getEmFolder()
        chimera_1_10_1_command = [('./scipion_installer',
                            '%s/chimera-1.13.1/bin/chimera' % SW_CH)]

        env.addPackage('chimera', version='1.13.1',
                       tar='chimera-1.13.1-linux_x86_64.tgz',
                       commands=chimera_1_10_1_command,
                       default=True)

        maxit_commands = [('make binary -j %d' % env.getProcessors(),
                            ['bin/maxit'])]

        env.addPackage('maxit', version='10.1',
                       tar=MAXIT_TAR,
                       url=MAXIT_URL,
                       commands=maxit_commands,
                       default=True)
pyworkflow.em.Domain.registerPlugin(__name__)
