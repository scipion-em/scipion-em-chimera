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

from ..protocols import ChimeraModelFromTemplate
from pwem.protocols.protocol_import import ProtImportPdb
from pwem.protocols.protocol_import.sequence import \
    ProtImportSequence
from pyworkflow.tests import *
from pwem.convert.sequence import indexToAlphabet


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
    UNIPROTID1 = "Q15116" # Programmed cell death protein 1
    UNIPROTID2 = "P01832" # Polymeric immunoglobulin receptor
    UNIPROTID3 = "Q2LC89"
    UNIPROTID4 = "Q9BQ51" # Programmed cell death 1 ligand 2
    UNIPROTID5 = "Q9NZQ7" # Programmed cell death 1 ligand 1

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
        protImportPDB.setObjLabel('import structure\n 1g03')
        self.launchProtocol(protImportPDB)
        structure3_PDB = protImportPDB.outputPdb
        return structure3_PDB

    def _importStructurePDBWoVol4(self):
        """Import the structure 3bp5 (mouse PD-1 and PD-L2 complex)"""
        args = {'inputPdbData': ProtImportPdb.IMPORT_FROM_ID,
                'pdbId': '3BP5'
                }
        protImportPDB = self.newProtocol(ProtImportPdb, **args)
        protImportPDB.setObjLabel('import structure\n 3bp5')
        self.launchProtocol(protImportPDB)
        structure4_PDB = protImportPDB.outputPdb
        return structure4_PDB

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

    def _importAminoacidSequence4(self):
        """Import the sequence derived from UniProtKB ID Q15116
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID1
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\nQ15116\n'
                                 'from UniProtKB')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("Q15116", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("Programmed cell death protein 1",
                         sequence.getDescription())
        self.assertEqual("MQIPQAPWPV", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence

    def _importAminoacidSequence5(self):
        """
         Import the sequence of chain A of atomic structure 3rrq.cif
         """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_STRUCTURE,
                'inputStructureSequence':
                    ProtImportSequence.IMPORT_STRUCTURE_FROM_ID,
                'pdbId': '3rrq',
                'inputStructureChain':
                    '{"model": 0, "chain": "A", "residues": 108}'
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\n from 3rrq\n'
                           'atomic structure')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence

        self.assertEqual("3rrq__0_A", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("", sequence.getDescription())
        self.assertEqual("NPPTFSPALL", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters)
        return sequence

    def _importAminoacidSequence6(self):
        """Import the sequence derived from UniProtKB ID P01832
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID2
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\nP01832\n'
                                 'from UniProtKB')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("P01832", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("Polymeric immunoglobulin receptor",
                         sequence.getDescription())
        self.assertEqual("MALFLLTCLL", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence

    def _importAminoacidSequence7(self):
        """Import the sequence derived from UniProtKB ID Q15116
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID3
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\nQ12LC89\n'
                                 'from UniProtKB')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("Q2LC89", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("PD-1-ligand 2",
                         sequence.getDescription())
        self.assertEqual("LQLHQIAALF", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence

    def _importAminoacidSequence8(self):
        """Import the sequence derived from UniProtKB ID Q9BQ51
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID4
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\nQ9BQ51\n'
                                 'from UniProtKB')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("Q9BQ51", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("Programmed cell death 1 ligand 2",
                         sequence.getDescription())
        self.assertEqual("MIFLLLMLSL", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence

    def _importAminoacidSequence9(self):
        """Import the sequence derived from UniProtKB ID Q9NZQ7
        """
        args = {'inputSequenceName': self.NAME,
                'inputProteinSequence':
                    ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': self.UNIPROTID5
                }
        protSequence = self.newProtocol(ProtImportSequence, **args)
        protSequence.setObjLabel('import aminoacid seq,\nQ9NZQ7\n'
                                 'from UniProtKB')
        self.launchProtocol(protSequence)
        sequence = protSequence.outputSequence
        self.assertEqual("Q9NZQ7", sequence.getId())
        self.assertEqual("User_Name", sequence.getSeqName())
        self.assertEqual("Programmed cell death 1 ligand 1",
                         sequence.getDescription())
        self.assertEqual("MRIFAVFIFM", sequence.getSequence()[:10])
        self.assertEqual("ACDEFGHIKLMNPQRSTVWY",
                         indexToAlphabet(sequence.getIsAminoacids(),
                                         sequence.getAlphabet()).letters
                         )
        return sequence


class TestChimeraModellerSearch(TestImportData):
    CHAIN1 = '{"model": 0, "chain": "B", "residues": 146}'  # Protein
    CHAIN2 = '{"model": 0, "chain": "A", "residues": 98}'  # Protein
    CHAIN3 = '{"model": 7, "chain": "A", "residues": 134}'  # Protein
    CHAIN4 = '{"model": 0, "chain": "A", "residues": 114}'
    CHAIN5 = '{"model": 0, "chain": "B", "residues": 191}'
    message = """*******************************************
This test requires human intervention. By default is disabled
so automatic checking do not fail but should be executed
if you want to check chimera modeler search protocol.

You may enable the protocol but setting

    DISABLE_TEST = False
    
HOW TO RUN THE TEST:

In the particular case of test testImportSequence5
you should select a possible template for the target, 
for example 5wt9, write in the line command "open 5wt9"  

1) when Chimerax opens select (in main menu)
Tools -> Sequence -> Modeller Comparative
2) a new pop up windows appears called 'Modeller Comparative'
3) template(s) are the Sequence alignments
   (in the case of test testImportSequence5, only the sequence alignment
   aligned_1.fasta)
4) choose as target(s) = seq...mutated
(In the particular case of test testImportChainFromStructureAndSequence4
select as target sequences Q15116 and Q2LC89 for aligned_1.fasta and 
aligned_2.fasta alignments, respectively. Do not forget select the 
basic option "Make multichain model from multichain template")
(in the case of test testImportSequence5, only the target sequence Q15116
in the aligned_1.fasta)

5) you will need a key for modeller website
6) press OK
7) In main Chimerax window you will see (botton) the progression of the process,
e.g. 1 of 5 models generated, 2 of 5 models generated, and so on
8) In command line rename the #id of the model in you are interested, 
for example from #3.1 to #4 type:
rename #3.1 id #4, enter>
9) In command line type:  scipionwrite #4 prefix model_3.1_,enter>
10) Close Chimerax 
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
                'inputSequence1': sequence1
                }
        prot1 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot1.setObjLabel('1_structure chain seq,\n and seq from '
                          'user file\n')
        if self.DISABLE_TEST:
            print(self.message)
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
                'inputSequence1': sequence2
                }
        prot2 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot2.setObjLabel('2_structure chain seq,\n and seq from '
                          'user file\n')
        if self.DISABLE_TEST:
            print(self.message)
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
                'inputSequence1': sequence3
                }
        prot3 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot3.setObjLabel('3_structure chain seq,\n and seq from '
                          'user file\n')
        if self.DISABLE_TEST:
            print(self.message)
        else:
            self.launchProtocol(prot3)

    def testImportChainFromStructureAndSequence4(self):
        """
        Import the aminoacid chains A and B from the structure 3bp5,
        and two target sequences (Q15116 and Q12LC89)
        """
        structure4_PDB = self._importStructurePDBWoVol4()
        sequence4 = self._importAminoacidSequence4()
        sequence5 = self._importAminoacidSequence5()
        sequence6 = self._importAminoacidSequence6()
        sequence7 = self._importAminoacidSequence7()
        sequence8 = self._importAminoacidSequence8()
        sequence9 = self._importAminoacidSequence9()

        args = {'pdbFileToBeRefined': structure4_PDB,
                'inputStructureChain': self.CHAIN4,
                'inputSequence1': sequence4,
                'optionForAligning1': 1,
                'inputSequencesToAlign1': [sequence5, sequence6],
                'additionalTargetSequence': True,
                'selectStructureChain': self.CHAIN5,
                'inputSequence2': sequence7,
                'optionForAligning2': 1,
                'inputSequencesToAlign2': [sequence8, sequence9]
                }
        prot3 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot3.setObjLabel('4_1 target sequences,\n '
                          'fitting a complex\n'
                          'additional seqs')
        if self.DISABLE_TEST:
            print(self.message)
        else:
            self.launchProtocol(prot3)

        alignmentFile1 = prot3._getExtraPath("aligned_1.fasta")
        alignmentFile2 = prot3._getExtraPath("aligned_2.fasta")

        args = {'pdbFileToBeRefined': structure4_PDB,
                'inputStructureChain': self.CHAIN4,
                'inputSequence1': sequence4,
                'optionForAligning1': 2,
                'inputYourOwnSequenceAlignment1': alignmentFile1,
                'additionalTargetSequence': True,
                'selectStructureChain': self.CHAIN5,
                'inputSequence2': sequence7,
                'optionForAligning2': 2,
                'inputYourOwnSequenceAlignment2': alignmentFile2,
                }
        prot4 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot4.setObjLabel('4_2 target sequences,\n '
                          'fitting a complex\n'
                          'alignment input')
        if self.DISABLE_TEST:
            print(self.message)
        else:
            self.launchProtocol(prot4)

    def testImportSequence5(self):
        """
        Import a target sequence Q15116 and look for templates
        """
        sequence4 = self._importAminoacidSequence4()
        args = {'addTemplate': False,
                'inputSequence1': sequence4,
                }
        prot3 = self.newProtocol(ChimeraModelFromTemplate, **args)
        prot3.setObjLabel('5 target sequence,\n '
                          'no template\n')
        if self.DISABLE_TEST:
            print(self.message)
        else:
            self.launchProtocol(prot3)