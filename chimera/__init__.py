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
import tempfile

from pathlib import Path

import pwem
import pyworkflow.utils as pwutils
from glob import glob

from pyworkflow import SPA, TOMO, MODELLING

from .constants import (CHIMERA_HOME, ALPHAFOLD_HOME, ALPHAFOLD_DATABASE_DIR,
                        chimeraTARs, V1_11_1, CHIMERA_FLATPAK_ID)
from .flatpak import is_installed, ask_install_scope
from pyworkflow.utils import redStr

__version__ = "4.0.1"
_logo = "chimerax_logo.png"
_references = ['Goddard2018']


class Plugin(pwem.Plugin):
    _homeVar = CHIMERA_HOME
    _pathVars = [CHIMERA_HOME]
    _supportedVersions = chimeraTARs.keys()
    _currentVersion = V1_11_1
    _fullVersion = 'chimerax-%s' % _currentVersion
    _processingField = [SPA, TOMO, MODELLING]
    answer = None

    def __init__(self):
        super().__init__()
        # Change package name to be chimerax
        # Ideally we could have changed the folder from chimera to chimerax --> chimera module is highly used by
        # other plugins and would require updating all of them.
        self._name = "chimerax"

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
        # return f"flatpak run {CHIMERA_FLATPAK_ID}"
        cmd = cls.getHome('bin', progName)
        return str(cmd)

    @classmethod
    def getPython(cls, progName="python*"):
        return cls.getProgram()
        # """ Return the program binary that will be used. """
        # path = glob(cls.getHome('bin', progName))
        # # todo only run this "vglrun in test mode
        # # return "vglrun " +  path[0]
        # return path[0]

    @classmethod
    def defineBinaries(cls, env):
        # env.showOnly = True
        from scipion.install.funcs import VOID_TGZ

        # cls.defineChimeraXInstallation(env, V1_1, default=True)
        cls.defineChimeraXInstallation(env,
                                       cls._currentVersion,
                                       default=True,
                                       tarDir=chimeraTARs[cls._currentVersion])
        print("ChimeraX installation defined.")
        # Scipion plugin for chimera.
        # It will depend on the version currently active
        pathToPlugin = os.path.join(os.path.dirname(__file__),
                                    "Bundles", "scipion")
        pathToBinary = cls.getProgram()
        print("Path to ChimeraX binary: %s" % pathToBinary)
        activeVersion = cls.getActiveVersion()
        installationFlagFile = "installed-%s" % activeVersion
        # flatpack can not access the /tmp directory, 
        # so we create the temporary file in the home directory.
        with tempfile.NamedTemporaryFile(dir=Path.home(),
                                         mode="w",
                                         delete=False,
                                         suffix=".cxc") as tmpFile:
            tmpFile.write("toolshed install QScore\n")
            tmpFile.write(f"devel install {pathToPlugin}\n")
            tmpFn = tmpFile.name

        installPluginsCommand = [(f"{pathToBinary} --nogui --exit {tmpFn} && "
                                  f"touch {installationFlagFile}",
                                  installationFlagFile)]
        import inspect

        env.addPackage('scipionchimera', version=__version__,
                       tar=VOID_TGZ,
                       default=True,
                       # needsProgs=["flatpak"],  # error message is not clear when flatpak is not installed, 
                       # so we check it in the installation function and exit with a clear message.
                       commands=installPluginsCommand)

    @classmethod
    def defineChimeraXInstallation(cls,
                                   env,
                                   version,
                                   default=False,
                                   tarDir=None):

        # check if flatpak is installed
        import shutil
        import sys


        if shutil.which("flatpak") is None:
            print(redStr("Flatpak is not installed"))
            print(redStr("Please install Flatpak (sudo apt install flatpak) "
                         "and try again."))
            sys.exit(1)

        # print("Flatpak is installed.")

        if is_installed(CHIMERA_FLATPAK_ID):
            print(
                f"{CHIMERA_FLATPAK_ID} already installed. No Binary installed")
            return

        if cls.answer is None:
            cls.answer = ask_install_scope()
        if cls.answer == "user":
            print("Installing ChimeraX for the current user.")
            flatpak_cmd = f"flatpak install --user ChimeraX-{version}.flatpak"
        elif cls.answer == "system":
            flatpak_cmd = f"sudo flatpak install -y ChimeraX-{version}.flatpak"
        else:
            print("Installation cancelled.")
            sys.exit(0)

        from scipion.install.funcs import \
            VOID_TGZ  # Local import to avoid having scipion-app installed when building the package.

        getchimera_script = os.path.join(
            os.path.dirname(__file__), "getchimera.py")

        extractionDir = finalDir = os.path.join("bin", "ChimeraX")
        if tarDir:
            extractionDir = os.path.join("..", tarDir, extractionDir)

        chimera_cmds = [
            ("pip install https://github.com/scipion-em/tk_html_widgets/archive/master.zip", []),
            ("""cd .. &&\
                python %s %s""" % (getchimera_script, version),
                "../ChimeraX-%s.flatpak" % version),
            (f"""cd .. &&\
                 {flatpak_cmd} &&\
                 mkdir -p chimerax-{version}/bin &&\
                 echo 'flatpak run {CHIMERA_FLATPAK_ID} $*'> chimerax-{version}/bin/ChimeraX &&\
                 chmod +x chimerax-{version}/bin/ChimeraX""", extractionDir)]
        print("ChimeraX installation commands: %s" % chimera_cmds)
        env.addPackage('chimerax', version=version,
                       tar=VOID_TGZ,
                       default=default,
                       commands=chimera_cmds,
                       )
