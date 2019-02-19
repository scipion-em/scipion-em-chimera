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


from os.path import exists
from collections import Counter
from pyworkflow.em.constants import (SYM_I222, SYM_I222r, SYM_In25, SYM_In25r,
                                     SYM_CYCLIC, SYM_DIHEDRAL, SYM_TETRAHEDRAL,
                                     SYM_OCTAHEDRAL, SCIPION_SYM_NAME)

from chimera.protocols import ProtContacts
from pyworkflow.tests import BaseTest, setupTestProject
from pyworkflow.em.protocol.protocol_import import ProtImportPdb

class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

class TestChimeraContact(TestImportBase):
    # protocol to test the chimera computed contacts between pairs
    # of chains

    def _importStructurePDB(self, pdbID='6b1t'):
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': pdbID
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import pdb\n %s' % pdbID)
        self.launchProtocol(protImportPDB)
        return protImportPDB.outputPdb


    def testContacts(self):
        #import sequence
        pdb1 = self._importStructurePDB('6b1t') # A
        args = {'pdbFileToBeRefined': pdb1,
                'symmetryGroup': SYM_I222,
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
                }

        protContacts = self.newProtocol(ProtContacts, **args)
        protContacts.setObjLabel('chain pairs')
        self.launchProtocol(protContacts)
        # TODO: more test are needed may be with an smaller sample
        c, conn = protContacts.prepareDataBase(drop=False)
        tableName = protContacts.getTableName()
        sqlCommand = """SELECT count(*) FROM {tableName}""".format(tableName=tableName)
        c.execute(sqlCommand)
        row = c.fetchone()
        self.assertEqual(int(row[0]), 488698)