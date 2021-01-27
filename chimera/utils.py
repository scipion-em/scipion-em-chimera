# Utils module for shared methods
import os
import sys

from pwem.viewers.viewer_chimera import chimeraPdbTemplateFileName, chimeraMapTemplateFileName, sessionFile

from chimera.Bundles.scipion.src.constants import *


def getEnvDictionary(prot):
    """ Returns a dictionary to pass environment variables to "communicate" with scipionchimera Bundle"""
    # Get the values
    _chimeraPdbTemplateFileName = os.path.abspath(prot._getExtraPath(chimeraPdbTemplateFileName))
    _chimeraMapTemplateFileName = os.path.abspath(prot._getExtraPath(chimeraMapTemplateFileName))
    _sessionFile = os.path.abspath(prot._getExtraPath(sessionFile))
    protId = prot.getObjId()

    # Populate the dictionary and return
    envDict = {CHIMERA_PDB_TEMPLATE_FILE_NAME: _chimeraPdbTemplateFileName % protId,
               CHIMERA_MAP_TEMPLATE_FILE_NAME: _chimeraMapTemplateFileName % protId,
               SESSIONFILE: _sessionFile,
               PROTID: str(prot.getObjId()),
               SCIPIONPYTHON: sys.executable}
    return envDict
