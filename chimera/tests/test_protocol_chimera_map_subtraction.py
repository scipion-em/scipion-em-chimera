# ***************************************************************************
# * Authors:    Marta Martinez (mmmtnez@cnb.csic.es)
# *             Roberto Marabini (roberto@cnb.csic.es)
# *
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
# ***************************************************************************/

# TODO: Fill in the next line
# protocol to test a chimerax map subtraction starting from two maps or
# a map and an additional map derived from an atomic structure.

from ..protocols import ChimeraProtOperate
from ..protocols import ChimeraSubtractionMaps
from pwem.protocols.protocol_import import (ProtImportPdb,
                                            ProtImportVolumes)
from ..constants import CHIMERA_CYCLIC
import re
# from pwem import Domain


from pyworkflow.tests import *
import os.path


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import map volumes and atomic structures(PDBx/mmCIF files)
    """
    pdbID = "5ni1"  # Haemoglobin atomic structure
    chainID = '{"model": 0, "chain": "A", "residues": 141}'
    removeResidues = '{"index": "6-10", "residues": "EEKSA"}'
    def _importVolume(self):
        args = {'filesPath': self.dsModBuild.getFile('volumes/emd_3488.map'),
                'samplingRate': 1.05,
                'setOrigCoord': False
                }
        protImportVol = self.newProtocol(ProtImportVolumes, **args)
        protImportVol.setObjLabel('import volume haemoglobin\n with default '
                                  'origin\n')
        self.launchProtocol(protImportVol)
        volume = protImportVol.outputVolume
        return volume

    def _importAtomStruct(self):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': self.pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n 5ni1')
        self.launchProtocol(protImportPDB)
        structure = protImportPDB.outputPdb
        return structure


class TestChimeraSubtractMap(TestImportData):
    """ Test the chimera subtraction map protocol
    """

    def testChimeraSubtract1(self):
        """ This test checks the subtraction in Chimera of a
        model-derived map from an imported map """
        print("Run Chimera subtraction of a model-derived map "
              "from the imported volume\n")

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure = self._importAtomStruct()

        # create auxiliary CMD file for chimera operate (model fitting)
        extraCommands = ""
        extraCommands += "move -52.50,-52.50,-51.88 model #3 " \
                         "coord #2\n"
        extraCommands += "fitmap #3 in #2\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate,
                                       **args)
        protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera1)
        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        try:
            result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d.getFileName()"
                          % protChimera1.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))


        # TODO: These steps of protChimera2 can not be performed in protChimera1 because
        # when the map is saved keep an inappropriate origin
        extraCommands = ""
        # extraCommands += "molmap #2 2.1 gridSpacing 1.05 modelId 3\n"
        extraCommands += "molmap #2 2.1 gridSpacing 1.05 replace false\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_\n"
        extraCommands += "exit\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': result,
                }
        protChimera2 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera2.setObjLabel('chimera operate\n pdb\n save '
                                 'model-derived map')
        self.launchProtocol(protChimera2)
        try:
            result = eval("protChimera2.DONOTSAVESESSION_Map__3_%06d.getFileName()"
                          % protChimera2.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera2.DONOTSAVESESSION_Map__3_%06d"
                      % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 0,
                'inputVolume2': result,
                }
        protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera3.setObjLabel('chimera subtract\n map -\n'
                                 'model-derived map\n')
        self.launchProtocol(protChimera3)

        # Dynamically defined name of the variable because it does depend on
        # the protocol ID
        try:
            result = eval("protChimera3.difference_Map__8_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval("protChimera3.filtered_Map__9_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

    def testChimeraSubtract2(self):
        """ This test checks the subtraction in Chimera of a
        model with control mutations from an imported map """
        print("Run Chimera subtraction of a model with control "
              "mutations from an imported volume\n")

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure = self._importAtomStruct()

        # create auxiliary CMD file for chimera operate (model fitting)
        extraCommands = ""
        extraCommands += "move -52.50,-52.50,-51.88 model #3 " \
                         "coord #2\n"
        extraCommands += "fitmap #3 in #2\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera1)
        try:
            result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d.getFileName()"
                 % protChimera1.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                'selectChain': True,
                'selectStructureChain': self.chainID,
                'removeResidues': True,
                'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2,
                'rangeDist': 100,
                'selectAreaMap': True,
                'radius': 1
                }
        protChimera2 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera2.setObjLabel('chimera subtract\n chainA-sym-derived map\n '
                                  'removed residues\nzone')
        self.launchProtocol(protChimera2)
        try:
            result = eval("protChimera2.chain_A_Atom_struct__4_%06d_cif.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera2.zone_Map__6_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera2.molmap_chainA_Map__7_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera2.difference_Map__8_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera2.filtered_Map__9_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                'selectChain': True,
                'selectStructureChain': self.chainID,
                'removeResidues': True,
                'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2,
                'rangeDist': 100
                }
        protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera3.setObjLabel('chimera subtract\n chainA-sym-derived map\n '
                                 'removed residues\nno zone')
        self.launchProtocol(protChimera3)
        try:
            result = eval("protChimera3.chain_A_Atom_struct__4_%06d_cif.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval("protChimera3.sym_Atom_struct__5_%06d_cif.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera3.molmap_chainA_Map__7_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera3.difference_Map__8_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera3.filtered_Map__9_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                'selectChain': True,
                'selectStructureChain': self.chainID,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2,
                'rangeDist': 100,
                'selectAreaMap': True,
                'radius': 1
                }
        protChimera4 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera4.setObjLabel('chimera subtract\n chainA-sym-derived map\n')
        self.launchProtocol(protChimera4)
        try:
            result = eval("protChimera4.chain_A_Atom_struct__4_%06d_cif.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval("protChimera4.sym_Atom_struct__5_%06d_cif.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval("protChimera4.zone_Map__6_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera4.molmap_chainA_Map__7_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera4.difference_Map__8_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera4.filtered_Map__9_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                'selectChain': True,
                'selectStructureChain': self.chainID,
                'selectAreaMap': True,
                'radius': 1
                }
        protChimera5 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera5.setObjLabel('chimera subtract\n chainA-derived map\n')
        self.launchProtocol(protChimera5)
        try:
            result = eval("protChimera5.chain_A_Atom_struct__4_%06d_cif.getFileName()"
                 % protChimera5.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera5.zone_Map__6_%06d.getFileName()"
                 % protChimera5.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera5.molmap_chainA_Map__7_%06d.getFileName()"
                 % protChimera5.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera5.difference_Map__8_%06d.getFileName()"
                 % protChimera5.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera5.filtered_Map__9_%06d.getFileName()"
                 % protChimera5.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                'selectChain': True,
                'selectStructureChain': self.chainID,
                }
        protChimera6 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera6.setObjLabel('chimera subtract\n chainA-derived map\nno zone')
        self.launchProtocol(protChimera6)
        try:
            result = eval("protChimera6.chain_A_Atom_struct__4_%06d_cif.getFileName()"
                 % protChimera6.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera6.molmap_chainA_Map__7_%06d.getFileName()"
                 % protChimera6.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera6.difference_Map__8_%06d.getFileName()"
                 % protChimera6.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera6.filtered_Map__9_%06d.getFileName()"
                 % protChimera6.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

    def testChimeraSubtract3(self):
        """ This test checks the subtraction in Chimera of a
        model with control mutations from an imported map """
        print("Run Chimera subtraction of a model with control "
              "mutations from an imported volume\n")

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure = self._importAtomStruct()

        # create auxiliary CMD file for chimera operate (model fitting)
        extraCommands = ""
        extraCommands += "move -52.50,-52.50,-51.88 model #3 " \
                         "coord #2\n"
        extraCommands += "fitmap #3 in #2\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera1)
        try:
            result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d.getFileName()"
                 % protChimera1.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                'removeResidues': True,
                'inputStructureChain': self.chainID,
                'residuesToRemove': self.removeResidues,
                }
        protChimera2 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera2.setObjLabel('chimera subtract\n atomStruct-derived map\n '
                                  'removed residues\nno zone')
        self.launchProtocol(protChimera2)
        try:
            result = eval("protChimera2.mutated_Atom_struct__3_%06d_cif.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera2.molmap_Map__7_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera2.difference_Map__8_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera2.filtered_Map__9_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval("protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
                      % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                }
        protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera3.setObjLabel('chimera subtract\n atomStruct-derived map\n '
                                 'no zone')
        self.launchProtocol(protChimera3)
        try:
            result = eval("protChimera3.Atom_struct__3_%06d_cif.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera3.molmap_Map__7_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera3.difference_Map__8_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval("protChimera3.filtered_Map__9_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

    def testChimeraSubtract4(self):
        """ This test checks the subtraction in Chimera of a
        model with control mutations from an imported map """
        print("Run Chimera subtraction of a model with control "
              "mutations from an imported volume\n")

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure = self._importAtomStruct()

        # create auxiliary CMD file for chimera operate (model fitting)
        extraCommands = ""
        extraCommands += "move -52.50,-52.50,-51.88 model #3 " \
                         "coord #2\n"
        extraCommands += "fitmap #3 in #2\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_\n"
        extraCommands += "exit\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera1)
        try:
            result = eval(
                "protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d.getFileName()"
                 % protChimera1.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # TODO: These steps of protChimera2 can not be performed in protChimera1 because

        # second way of obtaining the starting atom structure (it generates an appropriate
        # cif by applying symmetry)
        extraCommands = ""
        extraCommands += "sel #2/A,B\n"
        extraCommands += "save /tmp/chainA_B.cif format mmcif models #2 relModel #1 selectedOnly true\n"
        extraCommands += "open /tmp/chainA_B.cif\n"
        extraCommands += "scipionwrite #3 " \
                         "prefix DONOTSAVESESSION_A_B_\n"
        extraCommands += "exit\n"

        result = eval(
            "protChimera1.DONOTSAVESESSION_Atom_struct__3_%06d"
            % protChimera1.getObjId())
        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': result
                }
        protChimera2 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera2.setObjLabel('chimera operate\n pdb\n save chain A_B')
        self.launchProtocol(protChimera2)
        try:
            result = eval(
                "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d.getFileName()"
                 % protChimera2.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval(
            "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d"
            % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                # 'selectChain': True,
                # 'selectStructureChain': self.chainID,
                # 'removeResidues': True,
                # 'inputStructureChain': self.chainID,
                # 'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2,
                # 'selectAreaMap': True,
                # 'radius': 1
                }
        protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera3.setObjLabel('chimera subtract\n chainAB-sym-derived map\n')
        self.launchProtocol(protChimera3)
        try:
            result = eval(
                "protChimera3.Atom_struct__3_%06d_cif.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera3.sym_Atom_struct__5_%06d_cif.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera3.molmap_Map__7_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera3.difference_Map__8_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera3.filtered_Map__9_%06d.getFileName()"
                 % protChimera3.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval(
            "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d"
            % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                # 'selectChain': True,
                # 'selectStructureChain': self.chainID,
                # 'removeResidues': True,
                # 'inputStructureChain': self.chainID,
                # 'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2,
                'selectAreaMap': True,
                'radius': 1
                }
        protChimera4 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera4.setObjLabel('chimera subtract\n chainA-sym-derived map\n'
                                 'zone')
        self.launchProtocol(protChimera4)
        try:
            result = eval(
                "protChimera4.Atom_struct__3_%06d_cif.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera4.sym_Atom_struct__5_%06d_cif.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera4.zone_Map__6_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera4.molmap_Map__7_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera4.difference_Map__8_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera4.filtered_Map__9_%06d.getFileName()"
                 % protChimera4.getObjId())
        except:
            self.assertTrue(False,  "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval(
            "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d"
            % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                # 'selectChain': True,
                # 'selectStructureChain': self.chainID,
                'removeResidues': True,
                'inputStructureChain': self.chainID,
                'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2,
                'selectAreaMap': True,
                'radius': 1
                }
        protChimera5 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera5.setObjLabel('chimera subtract\n chainA-sym-derived map\n'
                                 'remove residues\nzone')
        self.launchProtocol(protChimera5)
        try:
            result = eval(
                "protChimera5.Atom_struct__3_%06d_cif.getFileName()"
                % protChimera5.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera5.sym_Atom_struct__5_%06d_cif.getFileName()"
                % protChimera5.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera5.zone_Map__6_%06d.getFileName()"
                % protChimera5.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera5.molmap_Map__7_%06d.getFileName()"
                % protChimera5.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera5.difference_Map__8_%06d.getFileName()"
                % protChimera5.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera5.filtered_Map__9_%06d.getFileName()"
                % protChimera5.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval(
            "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d"
            % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                # 'selectChain': True,
                # 'selectStructureChain': self.chainID,
                'removeResidues': True,
                'inputStructureChain': self.chainID,
                'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2
                }
        protChimera6 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera6.setObjLabel('chimera subtract\n chainAB-sym-derived map\n'
                                 'remove residues\nno zone')
        self.launchProtocol(protChimera6)
        try:
            result = eval(
                "protChimera6.Atom_struct__3_%06d_cif.getFileName()"
                % protChimera6.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera6.mutated_Atom_struct__3_%06d_cif.getFileName()"
                % protChimera6.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera6.sym_Atom_struct__5_%06d_cif.getFileName()"
                % protChimera6.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera6.molmap_Map__7_%06d.getFileName()"
                % protChimera6.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera6.difference_Map__8_%06d.getFileName()"
                % protChimera6.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera6.filtered_Map__9_%06d.getFileName()"
                % protChimera6.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera map subtraction
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval(
            "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d"
            % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'level': 0.217,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                # 'selectChain': True,
                # 'selectStructureChain': self.chainID,
                'removeResidues': True,
                'inputStructureChain': self.chainID,
                'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2
                }
        protChimera7 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera7.setObjLabel('chimera subtract\n chainAB-sym-derived map\n'
                                 'remove residues\nlevel selected\nno zone')
        self.launchProtocol(protChimera7)
        try:
            result = eval(
                "protChimera7.Atom_struct__3_%06d_cif.getFileName()"
                % protChimera7.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera7.mutated_Atom_struct__3_%06d_cif.getFileName()"
                % protChimera7.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera7.sym_Atom_struct__5_%06d_cif.getFileName()"
                % protChimera7.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera7.molmap_Map__7_%06d.getFileName()"
                % protChimera7.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera7.difference_Map__8_%06d.getFileName()"
                % protChimera7.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera7.filtered_Map__9_%06d.getFileName()"
                % protChimera7.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera mask
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval(
            "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d"
            % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'subtractOrMask': 1,
                'level': 0.217,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                # 'selectChain': True,
                # 'selectStructureChain': self.chainID,
                'removeResidues': True,
                'inputStructureChain': self.chainID,
                'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2
                }
        protChimera8 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera8.setObjLabel('chimera mask\n chainAB-sym-derived map\n'
                                 'remove residues\nlevel 0.217\nno zone')
        self.launchProtocol(protChimera8)
        try:
            result = eval(
                "protChimera8.Atom_struct__3_%06d_cif.getFileName()"
                % protChimera8.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera8.mutated_Atom_struct__3_%06d_cif.getFileName()"
                % protChimera8.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera8.sym_Atom_struct__5_%06d_cif.getFileName()"
                % protChimera8.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera8.molmap_Map__7_%06d.getFileName()"
                % protChimera8.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera8.difference_Map__8_%06d.getFileName()"
                % protChimera8.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera8.filtered_Map__9_%06d.getFileName()"
                % protChimera8.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        # protocol chimera mask
        extraCommands = "run(session, 'select all')\n"
        extraCommands += "run(session, 'exit')\n"
        result = eval(
            "protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3_%06d"
            % protChimera2.getObjId())
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'subtractOrMask': 1,
                'level': 2.0,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': result,
                # 'selectChain': True,
                # 'selectStructureChain': self.chainID,
                'removeResidues': True,
                'inputStructureChain': self.chainID,
                'residuesToRemove': self.removeResidues,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2
                }
        protChimera9 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera9.setObjLabel('chimera mask\n chainAB-sym-derived map\n'
                                 'remove residues\nlevel 2.0\nno zone')
        self.launchProtocol(protChimera9)
        try:
            result = eval(
                "protChimera9.Atom_struct__3_%06d_cif.getFileName()"
                % protChimera9.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera9.mutated_Atom_struct__3_%06d_cif.getFileName()"
                % protChimera9.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera9.sym_Atom_struct__5_%06d_cif.getFileName()"
                % protChimera9.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

        try:
            result = eval(
                "protChimera9.molmap_Map__7_%06d.getFileName()"
                % protChimera9.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera9.difference_Map__8_%06d.getFileName()"
                % protChimera9.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))
        try:
            result = eval(
                "protChimera9.filtered_Map__9_%06d.getFileName()"
                % protChimera9.getObjId())
        except:
            self.assertTrue(False, "There was a problem with the alignment")

        self.assertTrue(os.path.exists(result))

