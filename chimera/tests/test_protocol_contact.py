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


import os
from ..constants import (CHIMERA_I222, CHIMERA_I2n3,
                         CHIMERA_CYCLIC,
                         CHIMERA_DIHEDRAL_X,
                         CHIMERA_TETRAHEDRAL,
                         CHIMERA_OCTAHEDRAL)

from ..protocols import ChimeraProtContacts
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols.protocol_import import ProtImportPdb


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')


class TestImportData(TestImportBase):
    """ Import atomic structures(PDBx/mmCIF files)
    """

    def _importStructureFromPDBId(self, pdbID):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n%s' % pdbID)
        self.launchProtocol(protImportPDB)
        self.assertTrue(protImportPDB.outputPdb.getFileName())
        return protImportPDB.outputPdb

    def _importStructureFromFile(self, fileName):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(fileName)
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n%s' %
                                  os.path.basename(fileName))
        self.launchProtocol(protImportPDB)
        self.assertTrue(protImportPDB.outputPdb.getFileName())
        return protImportPDB.outputPdb


class TestChimeraContact(TestImportData):
    # protocol to test the chimera computed contacts between pairs
    # of chains

    def testContactsAsymetryC2(self):
        # import PDB; whole hemoglobin macromolecule with HEM groups as
        # independent chains
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1_HEM.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "A002": "HEM_A", '
                                  '"B": "chainB", "B002": "HEM_B", '
                                  '"C": "chainC", "C002": "HEM_C", '
                                  '"D": "chainD", "D002": "HEM_D"}',
                'applySymmetry': False
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_HEM\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 368)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 50  # 1 HEM_A       A002   #1 chainA      A
        # 57  # 1 HEM_B       B002   #1 chainB      B
        # 50  # 1 HEM_C       C002   #1 chainC      C
        # 58  # 1 HEM_D       D002   #1 chainD      D
        # 46  # 1 chainA      A     #1 chainB      B
        # 6  # 1 chainA      A     #1 chainC      C
        # 26  # 1 chainA      A     #1 chainD      D
        # 25  # 1 chainB      B     #1 chainC      C
        # 3  # 1 chainB      B     #1 chainD      D
        # 47  # 1 chainC      C     #1 chainD      D

    def testContactsAsymetryC2_b(self):
        # import PDB; whole hemoglobin macromolecule with HEM groups as
        # independent chains
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1_HEM.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "A002": "HEM_A", '
                                  '"B": "chainB", "B002": "HEM_B", '
                                  '"C": "chainC", "C002": "HEM_C", '
                                  '"D": "chainD", "D002": "HEM_D"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 1
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_HEM\nerror b\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 368)

        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 50  # 1 HEM_A       A002   #1 chainA      A
        # 57  # 1 HEM_B       B002   #1 chainB      B
        # 50  # 1 HEM_C       C002   #1 chainC      C
        # 58  # 1 HEM_D       D002   #1 chainD      D
        # 46  # 1 chainA      A     #1 chainB      B
        # 6  # 1 chainA      A     #1 chainC      C
        # 26  # 1 chainA      A     #1 chainD      D
        # 25  # 1 chainB      B     #1 chainC      C
        # 3  # 1 chainB      B     #1 chainD      D
        # 47  # 1 chainC      C     #1 chainD      D

    def testContactsAsymetryC2_c(self):
        # import PDB; whole hemoglobin macromolecule with HEM groups as
        # independent chains
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1_HEM.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "A002": "HEM_A", '
                                  '"B": "chainB", "B002": "HEM_B", '
                                  '"C": "chainC", "C002": "HEM_C", '
                                  '"D": "chainD", "D002": "HEM_D"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 0
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_HEM\nerror c\nno sym\ncontacts')
        try:
            self.launchProtocol(protContacts)
        except Exception as e:
            self.assertTrue(True)
            print("This test should return a error message as '"
                  " Error: Symmetry Order should be a positive integer.\n")

            return
        self.assertTrue(False)

    def testContactsSymC2_a(self):
        # import PDB; unit cell of hemoglobin macromolecule with HEM groups as
        # independent chains; sym C2
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1_unit_cell_HEM.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "A002": "HEM_A", '
                                  '"B": "chainB", "B002": "HEM_B"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_unit_cell_HEM\nsym C2\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 227)

        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 50  # 1.1 HEM_A       A002 #1.1 chainA      A
        # 57  # 1.1 HEM_B       B002 #1.1 chainB      B
        # 46  # 1.1 chainA      A   #1.1 chainB      B
        # 6  # 1.1 chainA      A   #1.2 chainA      A
        # 33  # 1.1 chainA      A   #1.2 chainB      B
        # 2  # 1.1 chainB      B   #1.2 chainB      B

    def testContactsSymC2_b(self):
        # import PDB; unit cell of hemoglobin macromolecule with HEM groups as
        # independent chains; sym C2
        # origin of coordinates different from center of symmetry
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "B": "chainB", '
                                  '"C": "chainC", "D": "chainD"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_CYCLIC,
                'symmetryOrder': 2
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1\nsym C2\ncenter of sym wrong\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 159)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 48  # 1 chainA      A     #1 chainB      B
        # 6  # 1 chainA      A     #1 chainC      C
        # 26  # 1 chainA      A     #1 chainD      D
        # 25  # 1 chainB      B     #1 chainC      C
        # 3  # 1 chainB      B     #1 chainD      D
        # 51  # 1 chainC      C     #1 chainD      D

    def testContactsAsymetryD4(self):
        # import PDB; whole molecule of thermosome from T. acidophilum (1a6d)
        # 8 independent chains, sym D4
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/1a6d_whole.pdb')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "up", "B": "up", "C": "up", "D": "up", '
                                  '"E": "up", "F": "down", "G": "down", "H": "down", '
                                  '"I": "down", "J": "up", "K": "up", "L": "up", '
                                  '"M": "down", "N": "down", "O": "down", "P": "down"}',
                'applySymmetry': False,
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1a6d_whole\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 408)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 56  # 1 down        F     #1 up          D
        # 56  # 1 down        G     #1 up          C
        # 56  # 1 down        H     #1 up          A
        # 56  # 1 down        I     #1 up          E
        # 46  # 1 down        M     #1 up          L
        # 46  # 1 down        N     #1 up          K
        # 46  # 1 down        O     #1 up          J
        # 46  # 1 down        P     #1 up          B

    def testContactsAsymetryD4_b(self):
        # import PDB; whole molecule of thermosome from T. acidophilum (1a6d)
        # 8 independent chains, sym D4
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/1a6d_whole.pdb')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "up", "B": "up", "C": "up", "D": "up", '
                                  '"E": "up", "F": "down", "G": "down", "H": "down", '
                                  '"I": "down", "J": "up", "K": "up", "L": "up", '
                                  '"M": "down", "N": "down", "O": "down", "P": "down"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_DIHEDRAL_X,
                'symmetryOrder': 1
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1a6d_whole\nerror b\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 408)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 56  # 1 down        F     #1 up          D
        # 56  # 1 down        G     #1 up          C
        # 56  # 1 down        H     #1 up          A
        # 56  # 1 down        I     #1 up          E
        # 46  # 1 down        M     #1 up          L
        # 46  # 1 down        N     #1 up          K
        # 46  # 1 down        O     #1 up          J
        # 46  # 1 down        P     #1 up          B

    def testContactsAsymetryD4_c(self):
        # import PDB; whole molecule of thermosome from T. acidophilum (1a6d)
        # 8 independent chains, sym D4
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/1a6d_whole.pdb')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "up", "B": "up", "C": "up", "D": "up", '
                                  '"E": "up", "F": "down", "G": "down", "H": "down", '
                                  '"I": "down", "J": "up", "K": "up", "L": "up", '
                                  '"M": "down", "N": "down", "O": "down", "P": "down"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_DIHEDRAL_X,
                'symmetryOrder': 0
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1a6d_whole\nerror c\nno sym\ncontacts')

        try:
            self.launchProtocol(protContacts)
        except Exception as e:
            self.assertTrue(True)
            print("This test should return a error message as '"
                  " Error: Symmetry Order should be a positive integer.\n")

            return
        self.assertTrue(False)

    def testContactsSymD4(self):
        # import PDB; unit cell of thermosome from T. acidophilum (1a6d)
        # 2 independent chains
        pdb1 = self._importStructureFromPDBId('1a6d')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "B": "chainB"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_DIHEDRAL_X,
                'symmetryOrder': 4
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1a6d_unit_cell\nsym D4\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 656)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 256  # 1.1 chainA      A   #1.1 chainB      B
        # 1  # 1.1 chainA      A   #1.2 chainA      A
        # 341  # 1.1 chainA      A   #1.2 chainB      B
        # 56  # 1.1 chainA      A   #1.4 chainA      A
        # 2  # 1.1 chainB      B   #1.2 chainB      B
        # 46  # 1.2 chainB      B   #1.4 chainB      B

    def testContactsAsymetryO(self):
        # import PDB; whole macromolecule of the cubic core of the pyruvate dehydrogenase
        # multienzyme complex (1eab); A. vinelandii; Ligands: Coenzyme A (CoA) and
        # 6,8-dimercapto-octanoic acid amide (LPM)
        # 24 independent chains
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/1eab_whole.pdb')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "up", "B": "up", "C": "down", "D": "down", '
                                  '"E": "down", "F": "down", "G": "down", "H": "down", '
                                  '"I": "down", "J": "down", "K": "up", "L": "up", '
                                  '"M": "up", "N": "up", "O": "up", "P": "up", '
                                  '"Q": "up", "R": "up", "S": "down", "T": "down", '
                                  '"U": "up", "V": "up", "W": "down", "X": "down"}',
                'applySymmetry': False,
                'symmetryGroup': CHIMERA_OCTAHEDRAL
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1eab_whole\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 240)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 60  # 1 down        C     #1 up          U
        # 60  # 1 down        I     #1 up          Q
        # 60  # 1 down        S     #1 up          A
        # 60  # 1 down        W     #1 up          K

    def testContactsSymO(self):
        # import PDB; unit cell of the cubic core of the pyruvate dehydrogenase
        # multienzyme complex (1eab); A. vinelandii; Ligands: Coenzyme A (CoA) and
        # 6,8-dimercapto-octanoic acid amide (LPM)
        # 1 independent chain
        pdb1 = self._importStructureFromPDBId('1eab')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_OCTAHEDRAL
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1eab_unit_cell\nsym O\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 274)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 133  # 1.1 chainA      A   #1.1 chainA      A
        # 274  # 1.1 chainA      A   #1.2 chainA      A
        # 272  # 1.2 chainA      A   #1.3 chainA      A

    def testContactsAsymetryT(self):
        # import PDB; tetrahedral oligomeric complex of GyrA N-terminal fragment
        # S. pneumoniae (6n1r);
        # 12 independent chains
        pdb1 = self._importStructureFromPDBId('6n1r')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "B": "chainB", "C": "chainC", '
                                  '"D": "chainD", "E": "chainE", "F": "chainF", '
                                  '"G": "chainG", "H": "chainH", "I": "chainI", '
                                  '"J": "chainJ", "K": "chainK", "L": "chainL"}',
                'applySymmetry': False,
                'symmetryGroup': CHIMERA_TETRAHEDRAL
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('6n1r\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 3641)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 186  # 1 chainA      A     #1 chainB      B
        # 210  # 1 chainA      A     #1 chainJ      J
        # 211  # 1 chainA      A     #1 chainL      L
        # 211  # 1 chainB      B     #1 chainI      I
        # 210  # 1 chainB      B     #1 chainK      K
        # 210  # 1 chainC      C     #1 chainF      F
        # 211  # 1 chainC      C     #1 chainH      H
        # 186  # 1 chainC      C     #1 chainI      I
        # 210  # 1 chainD      D     #1 chainE      E
        # 210  # 1 chainD      D     #1 chainG      G
        # 186  # 1 chainD      D     #1 chainJ      J
        # 186  # 1 chainE      E     #1 chainF      F
        # 211  # 1 chainE      E     #1 chainG      G
        # 210  # 1 chainF      F     #1 chainH      H
        # 186  # 1 chainG      G     #1 chainK      K
        # 186  # 1 chainH      H     #1 chainL      L
        # 211  # 1 chainI      I     #1 chainK      K
        # 210  # 1 chainJ      J     #1 chainL      L

    def testContactsSymT(self):
        # import PDB; unit cell of the tetrahedral aminopeptidase from P. horikoshii)
        # (1y0r); Ligands: Zn, As
        # 1 independent chain
        pdb1 = self._importStructureFromPDBId('1y0r')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA"}',
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_TETRAHEDRAL
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1y0r_unit_cell\nsym T\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 389)
        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 1116  # 1.1 chainA      A   #1.1 chainA      A
        # 256  # 1.1 chainA      A   #1.2 chainA      A
        # 133  # 1.1 chainA      A   #1.3 chainA      A
        # 133  # 1.3 chainA      A   #1.4 chainA      A

    def testContactsSymI222(self):
        """
        This test assesses contacts between any couple of proteins of the
        unit cell of an icosahedral virus and between any protein of the unit
        cell and any protein of the neighbour unit cells.
        """
        # import structure of the unit cell of a icosahedral virus
        pdb1 = self._importStructureFromPDBId('6b1t')  # A
        args = {'pdbFileToBeRefined': pdb1,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_I222,
                # symmetry group of the whole virus, once you have applied
                # symmetry to the unit cell
                'chainStructure': '{"A": "h1", "B": "h1", "C": "h1", '
                                  '"D": "h2", "E": "h2", "F": "h2", '
                                  '"G": "h3", "H": "h3", "I": "h3", '
                                  '"J": "h4", "K": "h4", "L": "h4", '
                                  '"M": "p", "N": "iiia", "O": "viiiO",'
                                  ' "P": "viiiP", "Q": "ix", '
                                  '"R": "ix", "S": "ix", '
                                  '"T": "ixb", "U": "vi", '
                                  '"V": "vi", "W": "vii", '
                                  '"X": "x", "Y": "vi"}'
                # labeling of unit cell chains
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('6b1t\nicosahedral virus\nsym I222\ncontacts ')
        self.launchProtocol(protContacts)
        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        # self.assertEqual(int(row[0]), 488698)
        self.assertEqual(int(row[0]), 19947)

        # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 54  # 1.1 h1          A   #1.1 h2          E
        # 8  # 1.1 h1          A   #1.1 h4          K
        # 2  # 1.1 h1          A   #1.1 iiia        N
        # 132  # 1.1 h1          A   #1.1 vi          Y
        # 32  # 1.1 h1          A   #1.1 viiiO       O
        # 4  # 1.1 h1          B   #1.1 iiia        N
        # 76  # 1.1 h1          B   #1.1 p           M
        # 1  # 1.1 h1          B   #1.1 vi          Y
        # 57  # 1.1 h1          B   #1.2 h1          B
        # 714  # 1.1 h1          B   #1.2 iiia        N
        # 4848  # 1.1 h1          B   #1.2 p           M
        # 83  # 1.1 h1          C   #1.1 h2          E
        # 100  # 1.1 h1          C   #1.1 h4          K
        # 10  # 1.1 h1          C   #1.1 iiia        N
        # 9  # 1.1 h1          C   #1.1 p           M
        # 90  # 1.1 h1          C   #1.1 vi          Y
        # 67  # 1.1 h1          C   #1.1 viiiO       O
        # 11  # 1.1 h1          C   #1.2 h1          C
        # 999  # 1.1 h1          C   #1.2 iiia        N
        # 1702  # 1.1 h1          C   #1.2 p           M
        # 8  # 1.1 h2          D   #1.1 h3          G
        # 87  # 1.1 h2          D   #1.1 ix          R
        # 3  # 1.1 h2          D   #1.1 vi          U
        # 114  # 1.1 h2          D   #1.1 vi          V
        # 60  # 1.1 h2          D   #1.1 vii         W
        # 27  # 1.1 h2          D   #1.1 viiiP       P
        # 31  # 1.1 h2          D   #1.1 x           X
        # 2  # 1.1 h2          E   #1.1 h4          J
        # 20  # 1.1 h2          E   #1.1 h4          K
        # 25  # 1.1 h2          E   #1.1 ix          Q
        # 13  # 1.1 h2          E   #1.1 ix          R
        # 13  # 1.1 h2          E   #1.1 ixb         T
        # 104  # 1.1 h2          E   #1.1 vi          U
        # 1  # 1.1 h2          E   #1.1 vi          V
        # 35  # 1.1 h2          E   #1.1 vii         W
        # 17  # 1.1 h2          E   #1.1 x           X
        # 90  # 1.1 h2          F   #1.1 h3          G
        # 33  # 1.1 h2          F   #1.1 h4          J
        # 19  # 1.1 h2          F   #1.1 h4          K
        # 34  # 1.1 h2          F   #1.1 ix          Q
        # 8  # 1.1 h2          F   #1.1 ixb         T
        # 104  # 1.1 h2          F   #1.1 vi          U
        # 107  # 1.1 h2          F   #1.1 vi          V
        # 70  # 1.1 h2          F   #1.1 viiiP       P
        # 51  # 1.1 h2          F   #1.1 x           X
        # 22  # 1.1 h3          G   #1.1 h4          J
        # 3  # 1.1 h3          G   #1.1 h4          K
        # 12  # 1.1 h3          G   #1.1 ix          Q
        # 2  # 1.1 h3          G   #1.1 ixb         T
        # 173  # 1.1 h3          G   #1.1 viiiP       P
        # 32  # 1.1 h3          H   #1.1 h4          J
        # 1  # 1.1 h3          H   #1.1 ix          Q
        # 127  # 1.1 h3          H   #1.1 ix          S
        # 40  # 1.1 h3          H   #1.1 viiiP       P
        # 81  # 1.1 h3          I   #1.1 ixb         T
        # 58  # 1.1 h4          J   #1.1 ix          Q
        # 8  # 1.1 h4          J   #1.1 ix          R
        # 3  # 1.1 h4          J   #1.1 ix          S
        # 133  # 1.1 h4          J   #1.1 viiiP       P
        # 10  # 1.1 h4          K   #1.1 iiia        N
        # 1  # 1.1 h4          K   #1.1 ix          S
        # 6  # 1.1 h4          K   #1.1 ixb         T
        # 162  # 1.1 h4          K   #1.1 viiiO       O
        # 86  # 1.1 h4          K   #1.1 viiiP       P
        # 32  # 1.1 h4          L   #1.1 iiia        N
        # 87  # 1.1 h4          L   #1.1 ix          Q
        # 36  # 1.1 h4          L   #1.1 viiiO       O
        # 3  # 1.1 h4          L   #1.1 viiiP       P
        # 60  # 1.1 iiia        N   #1.1 p           M
        # 1  # 1.1 iiia        N   #1.1 vi          Y
        # 84  # 1.1 iiia        N   #1.1 viiiO       O
        # 1331  # 1.1 iiia        N   #1.2 h1          A
        # 22  # 1.1 iiia        N   #1.2 h1          B
        # 776  # 1.1 iiia        N   #1.2 h1          C
        # 3900  # 1.1 iiia        N   #1.2 iiia        N
        # 424  # 1.1 iiia        N   #1.2 vi          Y
        # 886  # 1.1 iiia        N   #1.2 viiiO       O
        # 62  # 1.1 ix          Q   #1.1 ixb         T
        # 2  # 1.1 ix          R   #1.1 ixb         T
        # 72  # 1.1 ix          S   #1.1 ixb         T
        # 3  # 1.1 p           M   #1.2 h1          B
        # 18  # 1.1 p           M   #1.2 h1          C
        # 126  # 1.1 p           M   #1.2 iiia        N
        # 32  # 1.1 p           M   #1.2 p           M
        # 452  # 1.1 p           M   #1.2 viiiO       O
        # 4  # 1.1 vi          U   #1.1 vii         W
        # 13  # 1.1 vi          U   #1.1 x           X
        # 60  # 1.1 vi          V   #1.1 x           X
        # 9  # 1.1 vi          Y   #1.1 viiiO       O
        # 98  # 1.1 vi          Y   #1.2 iiia        N
        # 155  # 1.1 vi          Y   #1.2 p           M
        # 23  # 1.1 vii         W   #1.1 x           X
        # 29  # 1.1 viiiO       O   #1.2 h1          B
        # 84  # 1.1 viiiO       O   #1.2 h1          C
        # 28  # 1.1 viiiO       O   #1.2 p           M


    def testContactsSymI222_goodSym(self):
        """
        This test assesses contacts between any couple of proteins of the
        unit cell of an icosahedral virus and between any protein of the unit
        cell and any protein of the neighbour unit cells.
        """
        # import structure of the unit cell of a icosahedral virus
        pdb1 = self._importStructureFromPDBId('6b1t')  # A
        args = {'pdbFileToBeRefined': pdb1,
                'applySymmetry': True,
                'symmetryGroup': CHIMERA_I2n3,
                # symmetry group of the whole virus, once you have applied
                # symmetry to the unit cell
                'chainStructure': '{"A": "h1", "B": "h1", "C": "h1", '
                                  '"D": "h2", "E": "h2", "F": "h2", '
                                  '"G": "h3", "H": "h3", "I": "h3", '
                                  '"J": "h4", "K": "h4", "L": "h4", '
                                  '"M": "p", "N": "iiia", "O": "viiiO",'
                                  ' "P": "viiiP", "Q": "ix", '
                                  '"R": "ix", "S": "ix", '
                                  '"T": "ixb", "U": "vi", '
                                  '"V": "vi", "W": "vii", '
                                  '"X": "x", "Y": "vi"}'
                # labeling of unit cell chains
                }
        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('6b1t\nicosahedral virus\nsym I2n3\ncontacts ')
        self.launchProtocol(protContacts)
        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getView2Name()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertAlmostEqual(int(row[0]), 6968, delta=80)

        # # atoms, model_1, prot_1, chain_1,  model_2, prot_2, chain_2
        # 54  # 1.1 h1          A   #1.1 h2          E
        # 8  # 1.1 h1          A   #1.1 h4          K
        # 2  # 1.1 h1          A   #1.1 iiia        N
        # 132  # 1.1 h1          A   #1.1 vi          Y
        # 32  # 1.1 h1          A   #1.1 viiiO       O
        # 4  # 1.1 h1          B   #1.1 iiia        N
        # 76  # 1.1 h1          B   #1.1 p           M
        # 1  # 1.1 h1          B   #1.1 vi          Y
        # 3  # 1.1 h1          B   #1.4 h1          A
        # 68  # 1.1 h1          B   #1.4 h1          B
        # 1  # 1.1 h1          B   #1.4 p           M
        # 83  # 1.1 h1          C   #1.1 h2          E
        # 100  # 1.1 h1          C   #1.1 h4          K
        # 10  # 1.1 h1          C   #1.1 iiia        N
        # 9  # 1.1 h1          C   #1.1 p           M
        # 90  # 1.1 h1          C   #1.1 vi          Y
        # 67  # 1.1 h1          C   #1.1 viiiO       O
        # 41  # 1.1 h1          C   #1.4 h1          A
        # 33  # 1.1 h1          C   #1.4 h1          B
        # 8  # 1.1 h2          D   #1.1 h3          G
        # 87  # 1.1 h2          D   #1.1 ix          R
        # 3  # 1.1 h2          D   #1.1 vi          U
        # 114  # 1.1 h2          D   #1.1 vi          V
        # 60  # 1.1 h2          D   #1.1 vii         W
        # 27  # 1.1 h2          D   #1.1 viiiP       P
        # 31  # 1.1 h2          D   #1.1 x           X
        # 6  # 1.1 h2          D   #1.2 h2          D
        # 20  # 1.1 h2          D   #1.2 h2          F
        # 77  # 1.1 h2          D   #1.2 h3          G
        # 83  # 1.1 h2          D   #1.2 h3          I
        # 191  # 1.1 h2          D   #1.2 ixb         T
        # 2  # 1.1 h2          E   #1.1 h4          J
        # 20  # 1.1 h2          E   #1.1 h4          K
        # 25  # 1.1 h2          E   #1.1 ix          Q
        # 13  # 1.1 h2          E   #1.1 ix          R
        # 13  # 1.1 h2          E   #1.1 ixb         T
        # 104  # 1.1 h2          E   #1.1 vi          U
        # 1  # 1.1 h2          E   #1.1 vi          V
        # 35  # 1.1 h2          E   #1.1 vii         W
        # 17  # 1.1 h2          E   #1.1 x           X
        # 90  # 1.1 h2          F   #1.1 h3          G
        # 33  # 1.1 h2          F   #1.1 h4          J
        # 19  # 1.1 h2          F   #1.1 h4          K
        # 34  # 1.1 h2          F   #1.1 ix          Q
        # 8  # 1.1 h2          F   #1.1 ixb         T
        # 104  # 1.1 h2          F   #1.1 vi          U
        # 107  # 1.1 h2          F   #1.1 vi          V
        # 70  # 1.1 h2          F   #1.1 viiiP       P
        # 51  # 1.1 h2          F   #1.1 x           X
        # 20  # 1.1 h2          F   #1.2 h2          D
        # 9  # 1.1 h2          F   #1.2 h2          F
        # 22  # 1.1 h3          G   #1.1 h4          J
        # 3  # 1.1 h3          G   #1.1 h4          K
        # 12  # 1.1 h3          G   #1.1 ix          Q
        # 2  # 1.1 h3          G   #1.1 ixb         T
        # 173  # 1.1 h3          G   #1.1 viiiP       P
        # 32  # 1.1 h3          H   #1.1 h4          J
        # 1  # 1.1 h3          H   #1.1 ix          Q
        # 127  # 1.1 h3          H   #1.1 ix          S
        # 40  # 1.1 h3          H   #1.1 viiiP       P
        # 100  # 1.1 h3          H   #1.3 h3          H
        # 90  # 1.1 h3          H   #1.3 h3          I
        # 2  # 1.1 h3          H   #1.3 ix          S
        # 81  # 1.1 h3          I   #1.1 ixb         T
        # 5  # 1.1 h3          I   #1.2 ix          R
        # 58  # 1.1 h4          J   #1.1 ix          Q
        # 8  # 1.1 h4          J   #1.1 ix          R
        # 3  # 1.1 h4          J   #1.1 ix          S
        # 133  # 1.1 h4          J   #1.1 viiiP       P
        # 89  # 1.1 h4          J   #1.3 h3          I
        # 10  # 1.1 h4          K   #1.1 iiia        N
        # 1  # 1.1 h4          K   #1.1 ix          S
        # 6  # 1.1 h4          K   #1.1 ixb         T
        # 162  # 1.1 h4          K   #1.1 viiiO       O
        # 86  # 1.1 h4          K   #1.1 viiiP       P
        # 35  # 1.1 h4          K   #1.4 h1          A
        # 8  # 1.1 h4          K   #1.4 h1          B
        # 32  # 1.1 h4          L   #1.1 iiia        N
        # 87  # 1.1 h4          L   #1.1 ix          Q
        # 36  # 1.1 h4          L   #1.1 viiiO       O
        # 3  # 1.1 h4          L   #1.1 viiiP       P
        # 116  # 1.1 h4          L   #1.3 h3          I
        # 2  # 1.1 h4          L   #1.3 ixb         T
        # 54  # 1.1 h4          L   #1.4 h1          A
        # 88  # 1.1 h4          L   #1.4 h2          D
        # 92  # 1.1 h4          L   #1.4 h2          E
        # 183  # 1.1 h4          L   #1.4 ix          R
        # 60  # 1.1 iiia        N   #1.1 p           M
        # 1  # 1.1 iiia        N   #1.1 vi          Y
        # 84  # 1.1 iiia        N   #1.1 viiiO       O
        # 86  # 1.1 iiia        N   #1.4 h1          B
        # 96  # 1.1 iiia        N   #1.4 h1          C
        # 98  # 1.1 iiia        N   #1.4 iiia        N
        # 6  # 1.1 iiia        N   #1.4 p           M
        # 17  # 1.1 iiia        N   #1.4 vi          Y
        # 62  # 1.1 ix          Q   #1.1 ixb         T
        # 204  # 1.1 ix          Q   #1.3 h3          I
        # 44  # 1.1 ix          Q   #1.3 ixb         T
        # 4  # 1.1 ix          Q   #1.4 h2          D
        # 39  # 1.1 ix          Q   #1.4 ix          R
        # 2  # 1.1 ix          R   #1.1 ixb         T
        # 45  # 1.1 ix          R   #1.2 ixb         T
        # 72  # 1.1 ix          S   #1.1 ixb         T
        # 104  # 1.1 ix          S   #1.3 h3          H
        # 4  # 1.1 ix          S   #1.3 h3          I
        # 57  # 1.1 ix          S   #1.3 ix          S
        # 13  # 1.1 p           M   #1.4 h1          B
        # 12  # 1.1 p           M   #1.4 h1          C
        # 20  # 1.1 p           M   #1.4 iiia        N
        # 539  # 1.1 p           M   #1.4 p           M
        # 4  # 1.1 vi          U   #1.1 vii         W
        # 13  # 1.1 vi          U   #1.1 x           X
        # 60  # 1.1 vi          V   #1.1 x           X
        # 9  # 1.1 vi          Y   #1.1 viiiO       O
        # 23  # 1.1 vii         W   #1.1 x           X
        # 156  # 1.1 viiiO       O   #1.4 h1          A
        # 115  # 1.1 viiiO       O   #1.4 h1          B
        # 1  # 1.1 viiiO       O   #1.4 h1          C
        # 116  # 1.1 viiiO       O   #1.4 h2          E
        # 8  # 1.1 viiiO       O   #1.4 h2          F
        # 16  # 1.1 viiiP       P   #1.3 h3          G
        # 99  # 1.1 viiiP       P   #1.3 h3          I
