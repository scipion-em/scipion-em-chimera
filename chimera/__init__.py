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
from glob import glob
from .constants import (CHIMERA_HOME, ALPHAFOLD_HOME, ALPHAFOLD_DATABASE_DIR, 
                        V1_0, V1_1, V1_2_5, V1_3, V1_4, chimeraTARs)

__version__ = "3.3.4"
_logo = "chimerax_logo.png"
_references = ['Goddard2018']


class Plugin(pwem.Plugin):
    _homeVar = CHIMERA_HOME
    _pathVars = [CHIMERA_HOME]
    _supportedVersions = [V1_3]
    _currentVersion = V1_4  
    _fullVersion = 'chimerax-%s' % _currentVersion

    @classmethod
    def _defineVariables(cls):
        cls._defineEmVar(CHIMERA_HOME, cls._fullVersion)
        cls._defineVar(ALPHAFOLD_HOME, None)
        cls._defineVar(ALPHAFOLD_DATABASE_DIR, None)

    @classmethod
    def getEnviron(cls):
        environ = pwutils.Environ(os.environ)
        d = {}
        # d['PATH'] = cls.getHome('bin')
        d['PATH'] = cls.getHome('bin')
        if "REMOTE_MESA_LIB" in os.environ:
            d["LD_LIBRARY_PATH"] = os.environ['REMOTE_MESA_LIB']
        environ.update(d, position=pwutils.Environ.BEGIN)
        return environ

    @classmethod
    def runChimeraProgram(cls, program, args="", cwd=None, extraEnv=None):
        """ Internal shortcut function to launch chimera program. """
        env = cls.getEnviron()

        if extraEnv:
            env.update(extraEnv)

        pwutils.runJob(None, program, args, env=env, cwd=cwd)

    @classmethod
    def getProgram(cls, progName="ChimeraX"):
        """ Return the program binary that will be used. """
        cmd = cls.getHome('bin', progName)
        return str(cmd)

    @classmethod
    def getPython(cls, progName="python*"):
        """ Return the program binary that will be used. """
        path = glob(cls.getHome('bin', progName))
        # todo only run this "vglrun in test mode
        # return "vglrun " +  path[0]
        return path[0]

    @classmethod
    def isVersionActive(cls):
        return cls.getActiveVersion().startswith(V1_0)

    @classmethod
    def defineBinaries(cls, env):
        from scipion.install.funcs import VOID_TGZ

        #cls.defineChimeraXInstallation(env, V1_1, default=True)
        cls.defineChimeraXInstallation(env, cls._currentVersion, default=True, tarDir=chimeraTARs[cls._currentVersion])

        # Scipion plugin for chimera. It will depend on the version currently active
        pathToPlugin = os.path.join(os.path.dirname(__file__),
                                    "Bundles", "scipion")
        pathToBinary = cls.getProgram()

        activeVersion = cls.getActiveVersion()
        installationFlagFile = "installed-%s" % activeVersion
        ff = open("/tmp/kk.cxc", "w")
        ff.write(f'devel install {pathToPlugin}')
        ff.close()
        installPluginsCommand = [("%s --nogui --exit " \
                                  "/tmp/kk.cxc && touch %s" % (pathToBinary, installationFlagFile),
                                  [installationFlagFile])]

        env.addPackage('scipionchimera' , version='1.3',
                       tar=VOID_TGZ,
                       default=True,
                       commands=installPluginsCommand)

    @classmethod
    def defineChimeraXInstallation(cls, env, version, default=False, tarDir=None):
        from scipion.install.funcs import \
            VOID_TGZ  # Local import to avoid having scipion-app installed when building the package.

        getchimera_script = os.path.join(os.path.dirname(__file__),
                                        "getchimera.py")

        extractionDir = finalDir = os.path.join("bin", "ChimeraX")
        if tarDir:
            extractionDir = os.path.join("..", tarDir, extractionDir)

        chimera_cmds = [("cd .. && python %s %s" % (getchimera_script, version), "../ChimeraX-%s.tar.gz" %version),
                        ("cd .. && tar -xf ChimeraX-%s.tar.gz" % version, extractionDir)
                        ]

        if tarDir:
            chimera_cmds.append(("mv ../%s/* ." % tarDir,  finalDir))

        env.addPackage('chimerax', version=version,
                       tar=VOID_TGZ,
                       default=default,
                       commands=chimera_cmds,
                       )

