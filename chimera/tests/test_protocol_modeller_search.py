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


# protocol to test the chimera modeller search protocol starting from a
# pdb and a sequence. Here we are going to test the suitability of the
# protocol for uploading an alignment between the sequence of one of the
# structure chains (selected by the user) and the sequence that the user has
# introduced as input. Later on, the user has to send this alignment to
# Modeller server in order to get five independent models. The user has to
# select and save one of them.

from chimera.protocols import ChimeraModelFromTemplate
from pyworkflow.em.protocol.protocol_import import ProtImportPdb
from pyworkflow.em.protocol.protocol_import.sequence import \
    ProtImportSequence
from pyworkflow.tests import *
from pyworkflow.em.convert.sequence import indexToAlphabet


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')

class TestImportData(TestImportBase):
    """ Import atomic structures(PDBx/mmCIF files) and sequences
    """
    NAME = "User_Name"
    pdbID1 = "3lqd"  # Protein
    pdbID2 = "1aoi"  # DNA and protein


    def _importStructurePDBWoVol1(self):
        """Import the structure 3lqd (Protein structure)"""
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/3lqd.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n 3lqd')
        self.launchProtocol(protImportPDB)
        structure1_PDB = protImportPDB.outputPdb
        return structure1_PDB

    def _importStructurePDBWoVol2(self):
        """Import the structure 1aoi (Protein/DNA structure; nucleosome)"""
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1aoi.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n 1aoi')
        self.launchProtocol(protImportPDB)
        structure2_PDB = protImportPDB.outputPdb
        return structure2_PDB

    def _importStructurePDBWoVol3(self):
        """Import the structure 1g03 (NMR structure; 20 models)"""
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_FILES,
                'pdbFile': self.dsModBuild.getFile(
                    'PDBx_mmCIF/1g03.cif'),
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n 1aoi')
        self.launchProtocol(protImportPDB)
        structure3_PDB = protImportPDB.outputPdb
        return structure3_PDB

    def _importAminoacidSequence1(self):
        """Import the sequence derived from the mutation
        of 3lqd chain B sequence
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/3lqd_B_mutated.fasta')
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\n3lqd_B_mutated\n'
                                 'from file')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("Seq_3qld_B_mutated", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("Seq_3qld_B_mutated", sequence.getDescription())
        self.assertEqual("VHLSGEEKSA", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence

    def _importAminoacidSequence2(self):
        """Import the sequence derived from the mutation
        of 1aoi chain A sequence
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/1aoi_A_mutated.fasta')
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\n1aoi_A_mutated\n'
                                 'from file')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("Seq_1aoi_A_mutated", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("Seq_1aoi_A_mutated", sequence.getDescription())
        self.assertEqual("LATKAARKSS", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence

    def _importAminoacidSequence3(self):
        """Import the sequence derived from the mutation
        of 1g03 model 7 chain A sequence
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_FILES,
                'fileSequence': self.dsModBuild.getFile(
                    'Sequences/1g03__7_A_mutated.fasta')
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\n1g03__7_A_mutated\n'
                                 'from file')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("Seq_1g03__7_A_mutated", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("Seq_1g03__7_A_mutated", sequence.getDescription())
        self.assertEqual("PVMHPHGARR", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence

class TestChimeraModellerSearch(TestImportData):

    CHAIN1 = "[model: 0, chain: B, 148 residues]" # Protein
    CHAIN2 = "[model: 0, chain: A, 98 residues]" # Protein
    CHAIN3 = "[model: 7, chain: A, 134 residues]" # Protein
    message = """*******************************************
This test requires human intervention. By default is disabled
so automatic checking do not fail but should be executed
if you want to check chimera modeler search protocol.

You may enable the protocol but setting

    DISABLE_TEST = False
    
HOW TO RUN THE TEST:
1) when chimera opens select (in sequence window)
structure -> modeller (homology)
2) a new pop up windows appears called 'comparative model with modeller'
3) choose as target = seq...mutated
4) choose as template the only offered option 
5) you will need a key for modeller website
6) press OK
7) In main window you will see (botton) "porcentage of progress" wait untill is 100%
8) In command line type:  scipionwrite model #2.1 ,enter>
9) close chimera 
************************************************"""
    DISABLE_TEST = True 

    def testImportChainFromStructureAndSequence1(self):
        """
        Import the aminoacid chain B from the structure 3lqd and the sequence
        derived from the mutation of 3lqd chain B sequence
        """
        structure1_PDB = self._importStructurePDBWoVol1()
        sequence1 = self._importAminoacidSequence1()

        args = {'pdbFileToBeRefined': structure1_PDB,
                'inputStructureChain': self.CHAIN1,
                'inputSequence': sequence1
               }
        prot1 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot1.setObjLabel('1_structure chain seq,\n and seq from '
                          'user file\n')
        if self.DISABLE_TEST:
            print self.message
        else:
            self.launchProtocol(prot1)

    def testImportChainFromStructureAndSequence2(self):
        """
        Import the aminoacid chain A from the structure 1aoi (it contains
        both proteins and DNA) and the sequence
        derived from the mutation of 1aoi chain A sequence
        """
        structure2_PDB = self._importStructurePDBWoVol2()
        sequence2 = self._importAminoacidSequence2()

        args = {'pdbFileToBeRefined': structure2_PDB,
                'inputStructureChain': self.CHAIN2,
                'inputSequence': sequence2
               }
        prot2 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot2.setObjLabel('2_structure chain seq,\n and seq from '
                          'user file\n')
        if self.DISABLE_TEST:
            print self.message
        else:
            self.launchProtocol(prot2)

    def testImportChainFromStructureAndSequence3(self):
        """
        Import the aminoacid chain A from the structure 1g03, model 7, and the
        sequence derived from the mutation of 1g03, model 7, chain A sequence
        """
        structure3_PDB = self._importStructurePDBWoVol3()
        sequence3 = self._importAminoacidSequence3()

        args = {'pdbFileToBeRefined': structure3_PDB,
                'inputStructureChain': self.CHAIN3,
                'inputSequence': sequence3
                }
        prot3 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot3.setObjLabel('3_structure chain seq,\n and seq from '
                          'user file\n')
        if self.DISABLE_TEST:
            print self.message
        else:
            self.launchProtocol(prot3)
