# ***************************************************************************
# * Authors:   Roberto Marabini (roberto@cnb.csic.es)
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

# protocol to test the operation on PDB files
from os.path import exists
from collections import Counter

from ..protocols import ChimeraImportAtomStructAlphafold
from pyworkflow.tests import BaseTest, setupTestProject
import pwem.protocols as emprot


class TestBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

class TestChimeraAlphafoldImport(TestBase):

    def _auxiliaryFunctio(self, pdbID):
        pass 

    def testImportStructureFromEBIbyUniprotID(self):
        """Import atom struct prdiction by alphafold from EBI using uniprotid"""
        uniProtID = 'A0A087WSY6'
        args = {'source': ChimeraImportAtomStructAlphafold.IMPORT_FROM_EBI,
                'uniProtId': uniProtID
        }
        protImportAlpha = self.newProtocol(ChimeraImportAtomStructAlphafold, **args)
        protImportAlpha.setObjLabel('import from EBI with \n uniProtID %s' % uniProtID)
        self.launchProtocol(protImportAlpha)
        self.assertTrue(exists(protImportAlpha.A0A087WSY6.getFileName()))

    def testImportStructureFromEBIbyUniprotIDFail(self):
        """Import atom struct prdiction by alphafold from EBI using uniprotid"""
        uniProtID = 'AAAAAAAAA' # this ID does not exist
        args = {'source': ChimeraImportAtomStructAlphafold.IMPORT_FROM_EBI,
                'uniProtId': uniProtID
        }
        protImportAlpha = self.newProtocol(ChimeraImportAtomStructAlphafold, **args)
        protImportAlpha.setObjLabel('import from EBI with \n INVALID uniProtID %s' % uniProtID)
        with self.assertRaises(Exception) as context:
            self.launchProtocol(protImportAlpha)
        print("This test tests an exception. Therefore, the result should be 'OK'".upper())
        print("but you will see an ERROR on the screen".upper())

    def testImportStructureFromBlast(self):
        """Import atom struct predicted by alphafold
        using a sequence blast. saves the first returned model
        """
        # import sequence
        uniProtID = 'A0A087WSY6'
        args = {'inputSequenceName': 'seq_name',
                'inputProteinSequence': emprot.ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': uniProtID
                }

        prot = self.newProtocol(emprot.ProtImportSequence, **args)
        prot.setObjLabel('import aminoacid seq,\n from uniprot id')
        self.launchProtocol(prot)
        sequence = prot.outputSequence
        # save model2 since model1 in the axis
        extraCommands =  "run(session, 'alphafold fetch KVD15_HUMAN')\n"
        extraCommands += "run(session, 'scipionwrite #2 prefix DONOTSAVESESSION_')\n"
        extraCommands += "run(session, 'exit')\n"

        args = {'source': ChimeraImportAtomStructAlphafold.IMPORT_FROM_SEQ_BLAST,
                'inputSequence': sequence,
                'hideMessage': True,
                'extraCommands': extraCommands
                }
        prot1 = self.newProtocol(ChimeraImportAtomStructAlphafold, **args)
        prot1.setObjLabel('Alpha: Run blast \n')
        self.launchProtocol(prot1)
        self.assertTrue(exists(prot1.DONOTSAVESESSION_Atom_struct__2_000054.getFileName()))

    def testImportStructureFromColab(self):
        """Import atom struct predicted by alphafold
        using a sequence blast. saves the first returned model
        """
        message = """This test should be run interatively.
        I leave here this code so other developers know how to test it
        but remember you can NOT run the test. You need to open scipion3
        and creeate the protocols
        """
        INTERACTIVE = False
        if not INTERACTIVE:
            self.assertTrue(True)
            print(message)
            return
        # import sequence
        uniProtID = 'P0DKB3'
        args = {'inputSequenceName': 'seq_name',
                'inputProteinSequence': emprot.ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': uniProtID
                }

        prot = self.newProtocol(emprot.ProtImportSequence, **args)
        prot.setObjLabel('import from uniprotID')
        self.launchProtocol(prot)
        sequence = prot.outputSequence
        extraCommands = "run(session, 'exit')\n"

        args = {'source': ChimeraImportAtomStructAlphafold.IMPORT_REMOTE_ALPHAFOLD,
                'inputSequence': sequence,
                'hideMessage': True,
                'extraCommands': extraCommands
                }
        prot1 = self.newProtocol(ChimeraImportAtomStructAlphafold, **args)
        prot1.setObjLabel('execute colabs\n')
        self.launchProtocol(prot1)

    def testImportStructureLocal(self):
        message = """This test requires a valid alphafold instalation.
        So by default is not actived
        """
        INTERACTIVE = False
        if not INTERACTIVE:
            self.assertTrue(True)
            print(message)
            return

        # import sequence
        uniProtID = 'P0DKB3'
        args = {'inputSequenceName': 'seq_name',
                'inputProteinSequence': emprot.ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': uniProtID
                }

        prot = self.newProtocol(emprot.ProtImportSequence, **args)
        prot.setObjLabel('import from uniprotID')
        self.launchProtocol(prot)

        sequence = prot.outputSequence
        listSeq = [sequence, sequence]
        extraCommands = "run(session, 'exit')\n"
        args = {'source': ChimeraImportAtomStructAlphafold.IMPORT_LOCAL_ALPHAFOLD,
                'inputSequenceS': listSeq,
                'doGpu': True,
                'gpusToUse': 1,
                'extraFlags': '-c reduced_dbs',
                'hideMessage': True,
                'extraCommands': extraCommands
                }
        prot1 = self.newProtocol(ChimeraImportAtomStructAlphafold, **args)
        prot1.setObjLabel('execute local alphafold\n')
        self.launchProtocol(prot1)
