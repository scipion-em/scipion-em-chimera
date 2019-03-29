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
from collections import Counter
from pyworkflow.em.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL,
                                     SYM_OCTAHEDRAL, SCIPION_SYM_NAME)

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

    def testContactsAsymetry(self):
        # import PDB; whole hemoglobin macromolecule with HEM groups as
        # independent chains
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1_HEM.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "A002": "HEM_A", '
                                  '"B": "chainB", "B002": "HEM_B", '
                                  '"C": "chainC", "C002": "HEM_C", '
                                  '"D": "chainD", "D002": "HEM_D"}',
                'symmetryGroup': SYM_CYCLIC,
                'symmetryOrder': 1
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_HEM\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 736)

    def testContactsSymC2(self):
        # import PDB; unit cell of hemoglobin macromolecule with HEM groups as
        # independent chains; sym C2
        pdb1 = self._importStructureFromFile('PDBx_mmCIF/5ni1_unit_cell_HEM.cif')
        args = {'pdbFileToBeRefined': pdb1,
                'chainStructure': '{"A": "chainA", "A002": "HEM_A", '
                                  '"B": "chainB", "B002": "HEM_B"}',
                'symmetryGroup': SYM_CYCLIC,
                'symmetryOrder': 2
                }

        protContacts = self.newProtocol(ChimeraProtContacts, **args)
        protContacts.setObjLabel('5ni1_unit_cell_HEM\ncontacts')
        self.launchProtocol(protContacts)

        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 380)

    def testContactsSymI222(self):
        """
        This test assesses contacts between any couple of proteins of the
        unit cell of an icosahedral virus and between any protein of the unit
        cell and any protein of the neighbour unit cells.
        """
        # import structure of the unit cell of a icosahedral virus
        pdb1 = self._importStructureFromPDBId('6b1t') # A
        args = {'pdbFileToBeRefined': pdb1,
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
        # TODO: more test are needed may be with an smaller sample
        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        # self.assertEqual(int(row[0]), 488698)
        self.assertEqual(int(row[0]), 363616)