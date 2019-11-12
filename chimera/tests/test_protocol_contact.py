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
from pyworkflow.em.constants import (SYM_I222, SYM_I2n3,
                                     SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL,
                                     SYM_OCTAHEDRAL)

from chimera.protocols import ChimeraProtContacts
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pyworkflow.em.protocol.protocol_import import ProtImportPdb


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
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 736)

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
                'symmetryGroup': SYM_CYCLIC,
                'symmetryOrder': 1
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_HEM\nerror b\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 736)

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
                'symmetryGroup': SYM_CYCLIC,
                'symmetryOrder': 0
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_HEM\nerror c\nno sym\ncontacts')
        try:
            self.launchProtocol(protContacts)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " Error: Symmetry Order should be a positive integer.\n"

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
                'symmetryGroup': SYM_CYCLIC,
                'symmetryOrder': 2
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_unit_cell_HEM\nsym C2\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 380)

    def testContactsSymC2_b(self):
        # import PDB; unit cell of hemoglobin macromolecule with HEM groups as
        # independent chains; sym C2
        # origin of coordinates different from center of symmetry
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "B": "chainB", '
                                  '"C": "chainC", "D": "chainD"}',
                'applySymmetry': True,
                'symmetryGroup': SYM_CYCLIC,
                'symmetryOrder': 2
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1\nsym C2\ncenter of sym wrong\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 636)

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
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 816)

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
                'symmetryGroup': SYM_DIHEDRAL,
                'symmetryOrder': 1
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1a6d_whole\nerror b\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 816)

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
                'symmetryGroup': SYM_DIHEDRAL,
                'symmetryOrder': 0
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1a6d_whole\nerror c\nno sym\ncontacts')

        try:
            self.launchProtocol(protContacts)
        except Exception as e:
            self.assertTrue(True)
            print "This test should return a error message as '" \
                  " Error: Symmetry Order should be a positive integer.\n"

            return
        self.assertTrue(False)

    def testContactsSymD4(self):
        # import PDB; unit cell of thermosome from T. acidophilum (1a6d)
        # 2 independent chains
        pdb1 = self._importStructureFromPDBId('1a6d')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "B": "chainB"}',
                'applySymmetry': True,
                'symmetryGroup': SYM_DIHEDRAL,
                'symmetryOrder': 4
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1a6d_unit_cell\nsym D4\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 1302)

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
                'symmetryGroup': SYM_OCTAHEDRAL
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1eab_whole\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 480)

    def testContactsSymO(self):
        # import PDB; unit cell of the cubic core of the pyruvate dehydrogenase
        # multienzyme complex (1eab); A. vinelandii; Ligands: Coenzyme A (CoA) and
        # 6,8-dimercapto-octanoic acid amide (LPM)
        # 1 independent chain
        pdb1 = self._importStructureFromPDBId('1eab')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA"}',
                'applySymmetry': True,
                'symmetryGroup': SYM_OCTAHEDRAL
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1eab_unit_cell\nsym O\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 548)

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
                'symmetryGroup': SYM_TETRAHEDRAL
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('6n1r\nno sym\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 7282)


    def testContactsSymT(self):
        # import PDB; unit cell of the tetrahedral aminopeptidase from P. horikoshii)
        # (1y0r); Ligands: Zn, As
        # 1 independent chain
        pdb1 = self._importStructureFromPDBId('1y0r')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA"}',
                'applySymmetry': True,
                'symmetryGroup': SYM_TETRAHEDRAL
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('1y0r_unit_cell\nsym T\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 523)

    def testContactsSymI222(self):
        """
        This test assesses contacts between any couple of proteins of the
        unit cell of an icosahedral virus and between any protein of the unit
        cell and any protein of the neighbour unit cells.
        """
        # import structure of the unit cell of a icosahedral virus
        pdb1 = self._importStructureFromPDBId('6b1t') # A
        args = {'pdbFileToBeRefined': pdb1,
                'applySymmetry': True,
                'symmetryGroup': SYM_I222,
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
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        # self.assertEqual(int(row[0]), 488698)
        self.assertEqual(int(row[0]), 363096)

    def testContactsSymI222_goodSym(self):
        """
        This test assesses contacts between any couple of proteins of the
        unit cell of an icosahedral virus and between any protein of the unit
        cell and any protein of the neighbour unit cells.
        """
        # import structure of the unit cell of a icosahedral virus
        pdb1 = self._importStructureFromPDBId('6b1t') # A
        args = {'pdbFileToBeRefined': pdb1,
                'applySymmetry': True,
                'symmetryGroup': SYM_I2n3,
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
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertAlmostEqual(int(row[0]), 13081, delta=80)
