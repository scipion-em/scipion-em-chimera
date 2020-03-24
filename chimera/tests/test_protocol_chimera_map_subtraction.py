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
# protocol to test the chimera map subtraction starting from

from ..protocols import ChimeraProtOperate
from ..protocols import ChimeraSubtractionMaps
from pwem.protocols.protocol_import import (ProtImportPdb,
                                            ProtImportVolumes)
from ..constants import CHIMERA_CYCLIC
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
    FirstResidue = '{"residue": 6, "ASP"}'
    LastResidue = '{"residue": 10, "VAL"}'
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
        print("Run Chimera subtraction of a derived-model map "
              "from the imported volume\n")

        # Import Volume
        volume = self._importVolume()

        # import PDB
        structure = self._importAtomStruct()

        # create auxiliary CMD file for chimera operate (model fitting)
        extraCommands = ""
        extraCommands += "runCommand('move -52.50,-52.50,-51.88 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "prefix DONOTSAVESESSION_')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate,
                                       **args)
        protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera1)
        self.assertIsNotNone(
            protChimera1.DONOTSAVESESSION_Atom_struct__2.getFileName(),
                             "There was a problem with the alignment")
        self.assertTrue(os.path.exists(
            protChimera1.DONOTSAVESESSION_Atom_struct__2.getFileName()))


        # TODO: These steps of protChimera2 can not be performed in protChimera1 because
        # when the map is saved keep an inappropriate origin
        extraCommands = ""
        extraCommands += "runCommand('molmap #1 2.1 gridSpacing 1.05 modelId 2')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "prefix DONOTSAVESESSION_')\n"
        args = {'extraCommands': extraCommands,
                'pdbFileToBeRefined': protChimera1.DONOTSAVESESSION_Atom_struct__2
                }
        protChimera2 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera2.setObjLabel('chimera operate\n pdb\n save '
                                 'model-derived map')
        self.launchProtocol(protChimera2)
        self.assertIsNotNone(
            protChimera2.DONOTSAVESESSION_Map__2.getFileName(),
            "There was a problem with the alignment")
        self.assertTrue(os.path.exists(
            protChimera2.DONOTSAVESESSION_Map__2.getFileName()))


        # protocol chimera map subtraction
        extraCommands = "runCommand('select all')\n"
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 0,
                'inputVolume2': protChimera2.DONOTSAVESESSION_Map__2,
                }
        protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera3.setObjLabel('chimera subtract\n map -\n'
                                 'model-derived map\n')
        self.launchProtocol(protChimera3)
        self.assertIsNotNone(
            protChimera3.difference_Map__3.getFileName(),
            "There was a problem with the alignment")
        self.assertTrue(os.path.exists(
            protChimera3.filtered_Map__4.getFileName()))

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
        extraCommands += "runCommand('move -52.50,-52.50,-51.88 model #2 " \
                         "coord #1')\n"
        extraCommands += "runCommand('fitmap #2 #1')\n"
        extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
                         "prefix DONOTSAVESESSION_')\n"
        extraCommands += "runCommand('stop')\n"

        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'pdbFileToBeRefined': structure
                }
        protChimera1 = self.newProtocol(ChimeraProtOperate,
                                        **args)
        protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
                                 'fitted model')
        self.launchProtocol(protChimera1)
        self.assertIsNotNone(
            protChimera1.DONOTSAVESESSION_Atom_struct__2.getFileName(),
            "There was a problem with the alignment")
        self.assertTrue(os.path.exists(
            protChimera1.DONOTSAVESESSION_Atom_struct__2.getFileName()))

        # protocol chimera map subtraction
        extraCommands = "runCommand('select all')\n"
        args = {'extraCommands': extraCommands,
                'inputVolume': volume,
                'mapOrModel': 1,
                'resolution': 3.2,
                'pdbFileToBeRefined': protChimera1.DONOTSAVESESSION_Atom_struct__2,
                'selectChain': True,
                'selectStructureChain': self.chainID,
                'removeResidues': True,
                'firstResidueToRemove': self.FirstResidue,
                'lastResidueToRemove': self.LastResidue,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2,
                'rangeDist': 100,
                'selectAreaMap': True,
                'radius': 1
                }
        protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
        protChimera3.setObjLabel('chimera subtract\n chainA-sym-derived map\n '
                                  'removed residues')
        self.launchProtocol(protChimera3)
        self.assertIsNotNone(
            protChimera3.chain_A_Atom_struct__3.getFileName(),
            "There was a problem with the alignment")
        self.assertTrue(os.path.exists(
            protChimera3.zone_Map__6.getFileName()))
        self.assertTrue(os.path.exists(
            protChimera3.molmap_chainA_Map__7.getFileName()))
        self.assertTrue(os.path.exists(
            protChimera3.difference_Map__8.getFileName()))
        self.assertTrue(os.path.exists(
            protChimera3.filtered_Map__9.getFileName()))


    # def testChimeraSubtract3(self):
    #     """ This test checks the subtraction in Chimera of a
    #     model with control mutations from an imported map """
    #     print("Run Chimera subtraction of a model with control "
    #           "mutations from an imported volume\n")
    #
    #     # Import Volume
    #     volume = self._importVolume()
    #
    #     # import PDB
    #     structure = self._importAtomStruct()
    #
    #     # create auxiliary CMD file for chimera operate (model fitting)
    #     extraCommands = ""
    #     extraCommands += "runCommand('move -52.50,-52.50,-51.88 model #2 " \
    #                      "coord #1')\n"
    #     extraCommands += "runCommand('fitmap #2 #1')\n"
    #     extraCommands += "runCommand('scipionwrite model #2 refmodel #1 " \
    #                      "prefix DONOTSAVESESSION_')\n"
    #     extraCommands += "runCommand('stop')\n"
    #
    #     args = {'extraCommands': extraCommands,
    #             'inputVolume': volume,
    #             'pdbFileToBeRefined': structure
    #             }
    #     protChimera1 = self.newProtocol(ChimeraProtOperate,
    #                                     **args)
    #     protChimera1.setObjLabel('chimera operate\n volume and pdb\n save '
    #                              'fitted model')
    #     self.launchProtocol(protChimera1)
    #     self.assertIsNotNone(
    #         protChimera1.DONOTSAVESESSION_Atom_struct__2.getFileName(),
    #         "There was a problem with the alignment")
    #     self.assertTrue(os.path.exists(
    #         protChimera1.DONOTSAVESESSION_Atom_struct__2.getFileName()))
    #
    #     # TODO: These steps of protChimera2 can not be performed in protChimera1 because
    #     # when the map is saved keep an inappropriate origin
    #     extraCommands = ""
    #     extraCommands += "runCommand('split #1')\n"
    #     extraCommands += "runCommand('combine #1.1')\n"
    #     extraCommands += "runCommand('scipionwrite model #2 refmodel #0 " \
    #                      "prefix DONOTSAVESESSION_A_')\n"
    #     extraCommands += "runCommand('combine #1.1#1.2')\n"
    #     extraCommands += "runCommand('scipionwrite model #3 refmodel #0 " \
    #                      "prefix DONOTSAVESESSION_A_B_')\n"
    #     args = {'extraCommands': extraCommands,
    #             'pdbFileToBeRefined': protChimera1.DONOTSAVESESSION_Atom_struct__2
    #             }
    #     protChimera2 = self.newProtocol(ChimeraProtOperate,
    #                                     **args)
    #     protChimera2.setObjLabel('chimera operate\n pdb\n save '
    #                              'chain A and A_B')
    #     self.launchProtocol(protChimera2)
    #     self.assertIsNotNone(
    #         protChimera2.DONOTSAVESESSION_A_Atom_struct__2.getFileName(),
    #         "There was a problem with the alignment")
    #     self.assertTrue(os.path.exists(
    #         protChimera2.DONOTSAVESESSION_A_B_Atom_struct__3.getFileName()))
    #
    #     # protocol chimera map subtraction
    #     extraCommands = "runCommand('select all')\n"
    #     args = {'extraCommands': extraCommands,
    #             'inputVolume': volume,
    #             'mapOrModel': 1,
    #             'resolution': 3.2,
    #             'pdbFileToBeRefined': protChimera1.DONOTSAVESESSION_Atom_struct__2,
    #             'selectChain': True,
    #             'selectStructureChain': self.chainID,
    #             # 'removeResidues': True,
    #             # 'inputStructureChain': self.chainID,
    #             # 'firstResidueToRemove': self.FirstResidue,
    #             # 'lastResidueToRemove': self.LastResidue,
    #             # 'applySymmetry': True,
    #             # 'symmetryGroup': CHIMERA_CYCLIC,
    #             # 'symmetryOrder': 2,
    #             'selectAreaMap': True
    #             }
    #     protChimera3 = self.newProtocol(ChimeraSubtractionMaps, **args)
    #     protChimera3.setObjLabel('chimera subtract\n chainA-derived map\n '
    #                               'removed residues')
    #     self.launchProtocol(protChimera3)
    #     self.assertIsNotNone(
    #         protChimera3.chain_A_Atom_struct__3.getFileName(),
    #         "There was a problem with the alignment")
    #     self.assertTrue(os.path.exists(
    #         protChimera3.molmap_chainA_Map__4.getFileName()))
    #     self.assertTrue(os.path.exists(
    #         protChimera3.zone_Map__5.getFileName()))
    #     self.assertTrue(os.path.exists(
    #         protChimera3.difference_Map__6.getFileName()))
    #     self.assertTrue(os.path.exists(
    #         protChimera3.filtered_Map__7.getFileName()))
