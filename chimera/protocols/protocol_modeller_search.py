# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *
# **************************************************************************

import os
from pwem import *
from pwem.viewers.viewer_chimera import chimeraScriptFileName, Chimera, chimeraPdbTemplateFileName, \
    chimeraMapTemplateFileName, sessionFile
from pyworkflow import VERSION_1_2
from pyworkflow.utils import copyFile

from . import ChimeraProtBase

from pwem.convert.atom_struct import AtomicStructHandler

from pyworkflow.protocol.params import (PointerParam,
                                        StringParam,
                                        MultiPointerParam,
                                        BooleanParam,
                                        EnumParam,
                                        PathParam,
                                        FloatParam,
                                        IntParam)
from pwem.convert.sequence import (SequenceHandler,
                                   saveFileSequencesToAlign,
                                   alignClustalSequences,
                                   alignBioPairwise2Sequences,
                                   alignMuscleSequences)
from collections import OrderedDict

from ..constants import CLUSTALO, MUSCLE
from chimera import Plugin
from chimera.utils import getEnvDictionary


class ChimeraModelFromTemplate(ChimeraProtBase):
    """Protocol to model three-dimensional structures of proteins using Modeller.
        Execute command *scipionwrite #n [prefix stringAddedToFilename] from command line in order
        to transfer the selected
        pdb to scipion. Default value is model=#0,
        model refers to the pdb file."""
    _label = 'model from template'
    _program = ""
    _version = VERSION_1_2

    INFILE1 = "unaligned_1.fasta"
    OUTFILE1 = "aligned_1.fasta"
    INFILE2 = "unaligned_2.fasta"
    OUTFILE2 = "aligned_2.fasta"
    TWOSEQUENCES = 0
    MULTIPLESEQUENCES = 1
    ProgramToAlign1 = ['Bio.pairwise2', 'Clustal Omega', 'MUSCLE']
    ProgramToAlign2 = ['Clustal Omega', 'MUSCLE']
    OptionForAligning = ['None', 'Additional sequences to align',
                          'Provide your own sequence alignment']
    OptionForDataBase = ['PDB', 'NR']
    OptionForMatrix = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90',
                       'PAM30', 'PAM70', 'PAM250', 'IDENTITY']

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=False):
        form.addSection(label='Input')
        section = form.getSection('Input')

        section.addParam('addTemplate', BooleanParam,
                         default=True, label='Do you already have a template?',
                         help='"Yes": Option by default. Select this option in case '
                              'you already have a template to model your target '
                              'sequence.\n"No": Select this option if you want to '
                              'search for a template with which model your target '
                              'sequence. Generation of multimeric models is not '
                              'allowed selecting this option.\n')
        section.addParam('pdbFileToBeRefined', PointerParam,
                         pointerClass="AtomStruct", allowsNull=True,
                         important=True,
                         condition='addTemplate == True',
                         label='Atomic structure used as template',
                         help="PDBx/mmCIF file template used as basic atomic "
                              "structure to model your specific sequence.")
        section.addParam('inputStructureChain', StringParam,
                         label="Chain ", allowsNull=True, important=True,
                         condition='addTemplate == True',
                         help="Select a particular chain of the atomic "
                              "structure.")
        section.addParam('inputSequence1', PointerParam, pointerClass="Sequence",
                         label='Target sequence', allowsNull=True,
                         important=True,
                         help="Input the aminoacid sequence to align with the "
                              "structure template sequence.")
        section.addParam('dataBase', EnumParam,
                         choices=self.OptionForDataBase,
                         condition='addTemplate == False',
                         label="Protein sequence database:", default=0,
                         help="Select a protein sequence database to search "
                              "for templates:\nPDB: Experimentally determined structures "
                              "in the "
                              "Protein Data Bank.\nNR: NCBI 'non-redundant'database. "
                              "It contains GenBank translation proteins, PDB sequences, "
                              "SwissProt proteins + PIR + PRF. Since NR is much larger "
                              "than PDB, it takes longer to search.\n")
        section.addParam('similarityMatrix', EnumParam,
                         choices=self.OptionForMatrix,
                         condition='addTemplate == False',
                         label="Similarity matrix:", default=2,
                         help="Select a similarity matrix to use for alignment "
                              "scoring.\n")
        section.addParam('cutoffValue', FloatParam,
                         condition='addTemplate == False',
                         label="cutoff evalue:", default=1e-3,
                         help="Least significant expectation value needed to "
                              "qualify the retrieved element as a hit.\n")
        section.addParam('maxSeqs', IntParam,
                         condition='addTemplate == False',
                         label="Maximum number of sequences:", default=100,
                         help="Maximum number of sequences to retrieve "
                              "from the database.\n")
        section.addParam('optionForAligning1', EnumParam,
                         choices=self.OptionForAligning,
                         condition='addTemplate == True',
                         label="Options to improve the alignment:", default=0,
                         help="None: Option by default. Only the template and the "
                              "target sequences will be included in the alignment. "
                              "This option is recommendable when these two sequences "
                              "are very similar. Otherwise, select any of the two "
                              "additional options:\n"
                              "Additional sequences to align: Select this option "
                              "if you want to add some more sequences to accomplish "
                              "the alignment.\n"
                              "Provide your own sequence alignment: Your alignment"
                              "should include both the target and the template "
                              "sequences.\n")
        section.addParam('inputYourOwnSequenceAlignment1', PathParam,
                         pointerClass="File", allowsNull=False,
                         condition='addTemplate == True and optionForAligning1 == 2',
                         label='Sequence alignment input',
                         help="Input your own sequence alignment.\n"
                              "ChimeraX allowed formats accessible here: "
                              "https://www.cgl.ucsf.edu/chimerax/docs/user/commands/open.html#sequence ")
        section.addParam('inputSequencesToAlign1', MultiPointerParam,
                         pointerClass="Sequence", allowsNull=True,
                         condition='addTemplate == True and optionForAligning1 == 1',
                         label='Other sequences to align',
                         help="In case you need to load more sequences to "
                              "align, you can load them here.")
        section.addParam('inputProgramToAlign1_1', EnumParam,
                         choices=self.ProgramToAlign1,
                         label="Alignment tool for two sequences:", default=0,
                         condition='addTemplate == True and optionForAligning1 == 0',
                         help="Select a program to accomplish the sequence"
                              "alignment:\n\nBiophyton module "
                              "Bio.pairwise2 ("
                              "http://biopython.org/DIST/docs/api/"
                              "Bio.pairwise2-module.html). Built-in "
                              "program to align two "
                              "sequences. The global "
                              "alignment algorithm from the EMBOSS suite "
                              "has been implemented with match/mismatch "
                              "scores of 3/-1 and gap penalties "
                              "(open/extend) of "
                              "3/2.\n\nClustal Omega "
                              "program (http://www.clustal.org/omega/, "
                              "https://doi.org/10.1038/msb.2011.75): "
                              "Multiple sequence alignment tool. Install "
                              "clustalo if you choose this option for "
                              "the first time by 'sudo apt-get install "
                              "clustalo'.\n\nMUSCLE program stands for "
                              "MUltiple Sequence Comparison by "
                              "Log- Expectation("
                              "http://www.drive5.com/muscle/muscle.html, "
                              "https://doi.org/10.1093/nar/gkh340). "
                              "Install muscle if you choose this option "
                              "for the first time by 'sudo apt install "
                              "muscle'.")
        section.addParam('inputProgramToAlign2_1', EnumParam,
                         choices=self.ProgramToAlign2,
                         label="Multiple alignment tool:", default=0,
                         condition='addTemplate == True and optionForAligning1 == 1',
                         help="Select a program to accomplish the sequence"
                              "alignment:\n\nClustal Omega "
                              "program (http://www.clustal.org/omega/, "
                              "https://doi.org/10.1038/msb.2011.75): "
                              "Multiple sequence alignment tool. Install "
                              "clustalo if you choose this option for "
                              "the first time by 'sudo apt-get install "
                              "clustalo'.\n\nMUSCLE program stands for "
                              "MUltiple Sequence Comparison by "
                              "Log- Expectation("
                              "http://www.drive5.com/muscle/muscle.html, "
                              "https://doi.org/10.1093/nar/gkh340). "
                              "Install muscle if you choose this option "
                              "for the first time by 'sudo apt install "
                              "muscle'.")
        section.addParam('additionalTargetSequence', BooleanParam,
                         default=False,
                         condition='addTemplate == True',
                         label='Additional target sequence to include?',
                         help='Select YES if you want to add an additional '
                              'target sequence to model according a different '
                              'chain of the structure template. This '
                              'option is recommendable when you want to model '
                              'the two interacting elements of a particular complex'
                              ' at the same time.')
        section.addParam('selectStructureChain', StringParam,
                         condition='addTemplate == True and '
                                   'additionalTargetSequence == True',
                         label="Chain ", allowsNull=True, important=True,
                         help="Select a particular chain of the atomic "
                              "structure.")
        section.addParam('inputSequence2', PointerParam, pointerClass="Sequence",
                         condition='addTemplate == True and '
                                   'additionalTargetSequence == True',
                         label='Target sequence', allowsNull=True,
                         important=True,
                         help="Input the aminoacid sequence to align with the "
                              "structure template sequence.")
        section.addParam('optionForAligning2', EnumParam,
                         choices=self.OptionForAligning,
                         condition='addTemplate == True and '
                                   'additionalTargetSequence == True',
                         label="Options to improve the alignment:", default=0,
                         help="None: Option by default. Only the template and the "
                              "target sequences will be included in the alignment. "
                              "This option is recommendable when these two sequences "
                              "are very similar. Otherwise, select any of the two "
                              "additional options:\n"
                              "Additional sequences to align: Select this option "
                              "if you want to add some more sequences to accomplish "
                              "the alignment.\n"
                              "Provide your own sequence alignment: Your alignment"
                              "should include both the target and the template "
                              "sequences.\n")
        section.addParam('inputYourOwnSequenceAlignment2', PathParam,
                         pointerClass="File", allowsNull=False,
                         condition='addTemplate == True and '
                                   'optionForAligning2 == 2 and '
                                   'additionalTargetSequence == True',
                         label='Sequence alignment input',
                         help="Input your own sequence alignment.\n"
                              "ChimeraX allowed formats accessible here: "
                              "https://www.cgl.ucsf.edu/chimerax/docs/user/commands/open.html#sequence ")
        section.addParam('inputSequencesToAlign2', MultiPointerParam,
                         pointerClass="Sequence", allowsNull=True,
                         condition='addTemplate == True and '
                                   'optionForAligning2 == 1 and '
                                   'additionalTargetSequence == True',
                         label='Other sequences to align',
                         help="In case you need to load more sequences to "
                              "align, you can load them here.")
        section.addParam('inputProgramToAlign1_2', EnumParam,
                         choices=self.ProgramToAlign1,
                         label="Alignment tool for two sequences:", default=0,
                         condition='addTemplate == True and '
                                   'optionForAligning2 == 0 and '
                                   'additionalTargetSequence == True',
                         help="Select a program to accomplish the sequence"
                              "alignment:\n\nBiophyton module "
                              "Bio.pairwise2 ("
                              "http://biopython.org/DIST/docs/api/"
                              "Bio.pairwise2-module.html). Built-in "
                              "program to align two "
                              "sequences. The global "
                              "alignment algorithm from the EMBOSS suite "
                              "has been implemented with match/mismatch "
                              "scores of 3/-1 and gap penalties "
                              "(open/extend) of "
                              "3/2.\n\nClustal Omega "
                              "program (http://www.clustal.org/omega/, "
                              "https://doi.org/10.1038/msb.2011.75): "
                              "Multiple sequence alignment tool. Install "
                              "clustalo if you choose this option for "
                              "the first time by 'sudo apt-get install "
                              "clustalo'.\n\nMUSCLE program stands for "
                              "MUltiple Sequence Comparison by "
                              "Log- Expectation("
                              "http://www.drive5.com/muscle/muscle.html, "
                              "https://doi.org/10.1093/nar/gkh340). "
                              "Install muscle if you choose this option "
                              "for the first time by 'sudo apt install "
                              "muscle'.")
        section.addParam('inputProgramToAlign2_2', EnumParam,
                         choices=self.ProgramToAlign2,
                         label="Multiple alignment tool:", default=0,
                         condition='addTemplate == True and '
                                   'optionForAligning2 == 1 and '
                                   'additionalTargetSequence == True',
                         help="Select a program to accomplish the sequence"
                              "alignment:\n\nClustal Omega "
                              "program (http://www.clustal.org/omega/, "
                              "https://doi.org/10.1038/msb.2011.75): "
                              "Multiple sequence alignment tool. Install "
                              "clustalo if you choose this option for "
                              "the first time by 'sudo apt-get install "
                              "clustalo'.\n\nMUSCLE program stands for "
                              "MUltiple Sequence Comparison by "
                              "Log- Expectation("
                              "http://www.drive5.com/muscle/muscle.html, "
                              "https://doi.org/10.1093/nar/gkh340). "
                              "Install muscle if you choose this option "
                              "for the first time by 'sudo apt install "
                              "muscle'.")
        section.addParam('extraCommands', StringParam,
                          default = '',
                          condition = 'False',
                          label = 'Extra commands for chimera viewer',
                          help = "Add extra commands in cmd file. Use for testing")
        form.addSection(label='Help')
        form.addLine("Step 1:\nIn the sequence window your target "
                         "sequence (and other additional sequences that you "
                         "want to use in  the alignment) will appear aligned to "
                         "the template's sequence. Select in the sequence window "
                         "menu:\nTools -> Sequence -> Modeller Comparative;\nA new  "
                         "window for Comparative Modeling with Modeller will "
                         "appear. Select your specific template(s) as the Sequence "
                         "alignments and the target(s)sequence as the sequence "
                         "to be modeled"
                         + '''
        . To run Modeller via web service 
        write the Modeller license key supplied (Academic user can 
        register free of charge to receive a license key). Finally, press OK.
        \nWAITING TIME: (you may see the status of your job in chimera main 
        window, lower left corner.)\n\nStep 2:\nWhen the process finished, 
        5 models will 
        be automatically superimposed onto the template and model scores
        will appear in Modeller Results window. In Chimera Model panel 
        you will have: #1 (coordinate axes); #2 (
        template); #3.1 to 3.5 (models).Choose the one you like the best, 
        for example model #3.1. To save it in Scipion, we need to change the 
        model ID. In Chimera main menu: Favorites -> Command Line, write 
        *rename #3.1 id #4*. Then, you will see in Model panel 
        that selected model #3.1 renamed to #3. Save it 
        as first guess in Scipion by executing the Chimera command 
        *scipionwrite [model] #n [prefix XX]*. In our example 
        *scipionwrite #4 pefix model_3_1_*.\n 
        When you use the command line scipionwrite, the Chimera session will 
        be saved by default. Additionally, you can save the Chimera session 
        whenever you want by executing the command *scipionss*. You will be 
        able to restore the saved session by using the protocol chimera 
        restore session (SCIPION menu: Tools/Calculators/chimera restore 
        session). Once you have save your favorite model you can press 
        Quit in the Modeller Results window.''')

    # --------------------------- INSERT steps functions --------------------

    def prerequisitesStep(self):
        if self.addTemplate:

            # read PDB
            fileName = self._readPDB()

            # get pdb sequence
            import json
            chainIdDict = json.loads(self.inputStructureChain.get())

            userSeq = self.inputSequence1.get()  # SEQ object from Scipion

            inFile = self.INFILE1
            outFile = self.OUTFILE1

            addSeq = self.optionForAligning1.get()

            yourAlignment = self.inputYourOwnSequenceAlignment1.get()

            inputSeqAlign = self.inputSequencesToAlign1

            programToAlign1 = self.inputProgramToAlign1_1

            programToAlign2 = self.inputProgramToAlign2_1

            self.prePreRequisites(fileName, chainIdDict, userSeq,
                                  inFile, outFile, addSeq, yourAlignment,
                                  inputSeqAlign, programToAlign1,
                                  programToAlign2)

            self.selectedChain1 = self.selectedChain

            if self.additionalTargetSequence.get() is True:
                chainIdDict = json.loads(self.selectStructureChain.get())

                userSeq = self.inputSequence2.get()  # SEQ object from Scipion

                inFile = self.INFILE2
                outFile = self.OUTFILE2

                addSeq = self.optionForAligning2.get()

                yourAlignment = self.inputYourOwnSequenceAlignment2.get()

                inputSeqAlign = self.inputSequencesToAlign2

                programToAlign1 = self.inputProgramToAlign1_2

                programToAlign2 = self.inputProgramToAlign2_2

                self.prePreRequisites(fileName, chainIdDict, userSeq,
                                      inFile, outFile, addSeq, yourAlignment,
                                      inputSeqAlign, programToAlign1,
                                      programToAlign2)

            self.selectedChain2 = self.selectedChain

        else:
            userSeq = self.inputSequence1.get()  # SEQ object from Scipion
            # get target sequence imported by the user
            outFile = self.OUTFILE1
            self.targetSeqID1 = self.preTemplate(userSeq, outFile)

    def prePreRequisites(self, fileName, chainIdDict, userSeq, inFile,
                         outFile, addSeq, yourAlignment, inputSeqAlign,
                         programToAlign1, programToAlign2):
        # get sequence of structure chain with id chainId (selected by the user)
        self.selectedModel = chainIdDict['model']
        self.selectedChain = chainIdDict['chain']
        # self.selectedModel = chainId.split(',')[0].split(':')[1].strip()
        # self.selectedChain = chainId.split(',')[1].split(':')[1].strip()
        print("Selected chain: %s from model: %s from structure: %s" \
              % (self.selectedChain, self.selectedModel,
                 os.path.basename(fileName)))

        # Bio.Seq.Seq object
        structureSeq = self.structureHandler.getSequenceFromChain(
            self.selectedModel, self.selectedChain)

        # obtain a seqID for our PDB sequence
        structSeqID = self.structureHandler.getFullID(self.selectedModel,
                                                      self.selectedChain)
        # END PDB sequence

        # start user imported target sequence
        # get target sequence imported by the user

        targetSeqID = userSeq.getId()  # ID associated to SEQ object (str)
        userSequence = userSeq.getSequence()  # sequence associated to
        # that SEQ object (str)
        # transformation of this sequence (str) in a Bio.Seq.Seq object:
        seqHandler = SequenceHandler(userSequence,
                                     isAminoacid=userSeq.getIsAminoacids())
        targetSeq = seqHandler._sequence  # Bio.Seq.Seq object

        # creation of Dic of IDs and sequences
        SeqDic = OrderedDict()
        SeqDic[structSeqID] = structureSeq
        SeqDic[targetSeqID] = targetSeq

        # align sequences and save them to disk, -this will be chimera input-
        # get all sequences in a fasta file
        inFile = self._getInFastaSequencesFile(inFile)
        outFile = self._getOutFastaSequencesFile(outFile)

        # get the alignment of sequences
        if addSeq == 0:
            saveFileSequencesToAlign(SeqDic, inFile)
            inputSeqAlign = None
            if programToAlign1.get() == \
                    self.ProgramToAlign1.index('Bio.pairwise2'):
                # Only the two first sequences will be included in the alignment
                self.alignment = alignBioPairwise2Sequences(
                    structSeqID, structureSeq,
                    targetSeqID, targetSeq,
                    outFile)
            else:
                # All the sequences will be included in the alignment
                if programToAlign1.get() == \
                        self.ProgramToAlign1.index('Clustal Omega'):
                    cline = alignClustalSequences(inFile, outFile)
                else:
                    cline = alignMuscleSequences(inFile, outFile)
                args = ''
                self.runJob(cline, args)
        elif addSeq == 1:
            # if there are additional sequences imported by the user
            if inputSeqAlign is not None:
                for seq in inputSeqAlign:
                    seq = seq.get()
                    ID = seq.getId()
                    sequence = seq.getSequence()
                    seqHandler = SequenceHandler(sequence,
                                                 isAminoacid=seq.getIsAminoacids())
                    otherSeq = seqHandler._sequence  # Bio.Seq.Seq object
                    SeqDic[ID] = otherSeq

            # align sequences and save them to disk, -this will be chimera input-
            # get all sequences in a fasta file
            # inFile = self._getInFastaSequencesFile()
            saveFileSequencesToAlign(SeqDic, inFile)
            # outFile = self._getOutFastaSequencesFile()

            # All the sequences will be included in the alignment
            if programToAlign2 == self.ProgramToAlign2.index(
                    'Clustal Omega'):
                cline = alignClustalSequences(inFile, outFile)
            else:
                cline = alignMuscleSequences(inFile, outFile)
            args = ''
            self.runJob(cline, args)

        else:
            aligmentFile = os.path.basename(yourAlignment)
            outFile = os.path.join(self._getExtraPath(), aligmentFile)
            copyFile(yourAlignment, outFile)

    def preTemplate(self, userSeq, outFile):
        userSequence = userSeq.getSequence()  # sequence associated to
        # that SEQ object (str)
        targetSeqID = userSeq.getId()  # ID associated to SEQ object (str)
        # transformation of this sequence (str) in a Bio.Seq.Seq object:
        seqHandler = SequenceHandler(userSequence,
                                     isAminoacid=userSeq.getIsAminoacids())
        targetSeq = seqHandler._sequence  # Bio.Seq.Seq object
        # creation of Dic of IDs and sequences
        SeqDic = OrderedDict()
        SeqDic[targetSeqID] = targetSeq
        outFile = self._getOutFastaSequencesFile(outFile)
        saveFileSequencesToAlign(SeqDic, outFile)
        return targetSeqID

    def _readPDB(self):
        self.structureHandler = AtomicStructHandler()
        fileName = os.path.abspath(self.pdbFileToBeRefined.get(
        ).getFileName())
        self.structureHandler.read(fileName)
        return fileName

    def _getInFastaSequencesFile(self, inFile):
        INFILENAME = self._getTmpPath(inFile)
        return os.path.abspath(INFILENAME)

    def _getOutFastaSequencesFile(self, outFile):
        OUTFILENAME = self._getExtraPath(outFile)
        return os.path.abspath(OUTFILENAME)

    def runChimeraStep(self):
        # building script file including the coordinate axes and the input
        # volume with samplingRate and Origin information
        f = open(self._getTmpPath(chimeraScriptFileName), "w")
        # building coordinate axes

        dim = 150  # eventually we will create a PDB library that
        # computes PDB dim
        sampling = 1.

        tmpFileName = os.path.abspath(self._getTmpPath("axis_input.bild"))
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=tmpFileName,
                                         sampling=sampling)
        f.write("open %s\n" % tmpFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates

        # input vol with its origin coordinates
        pdbModelCounter = 1
        if (not self.addTemplate and
                self.inputSequence1.get() is not None and
                self._getOutFastaSequencesFile is not None):
            alignmentFile1 = self._getOutFastaSequencesFile(self.OUTFILE1)
            f.write("open %s\n" % alignmentFile1)
            f.write("blastprotein %s:%s database %s matrix %s "
                    "cutoff %.3f maxSeqs %d log true\n" %
                    (alignmentFile1.split("/")[-1], self.targetSeqID1,
                     self.OptionForDataBase[int(self.dataBase)],
                     self.OptionForMatrix[int(self.similarityMatrix)],
                     self.cutoffValue, self.maxSeqs))

        if (hasattr(self, 'pdbFileToBeRefined') and
                self.pdbFileToBeRefined.get() is not None):
            pdbModelCounter += 1
            pdbFileToBeRefined = self.pdbFileToBeRefined.get()
            f.write("open %s\n" % os.path.abspath(
                pdbFileToBeRefined.getFileName()))
            if pdbFileToBeRefined.hasOrigin():
                x, y, z = (pdbFileToBeRefined.getOrigin().getShifts())
                f.write("move %0.2f,%0.2f,%0.2f model #%d "
                        "coord #0\n" % (x, y, z, pdbModelCounter))

        # Alignment of sequence and structure
        if (hasattr(self, 'inputSequence1') and
                hasattr(self, 'inputStructureChain')):
            if (self.inputSequence1.get() is not None and
                self.inputStructureChain.get() is not None):
                pdbModelCounter = 2
                if str(self.selectedModel) != '0':
                    f.write("select #%s.%s/%s\n"
                            % (pdbModelCounter,
                               str(self.selectedModel + 1),
                               str(self.selectedChain1)))
                else:
                    f.write("select #%s/%s\n"
                            % (pdbModelCounter,
                               str(self.selectedChain1)))

                if self._getOutFastaSequencesFile is not None:
                    alignmentFile1 = self._getOutFastaSequencesFile(self.OUTFILE1)
                    f.write("open %s\n" % alignmentFile1)
                    f.write("sequence disassociate #%s %s\n" %
                            (pdbModelCounter,
                            alignmentFile1.split("/")[-1]))
                    if str(self.selectedModel) != '0':
                        f.write("sequence associate #%s.%s/%s %s:1\n" %
                                (pdbModelCounter,
                                 str(self.selectedModel + 1),
                                 str(self.selectedChain1),
                                 alignmentFile1.split("/")[-1]))
                    else:
                        f.write("sequence associate #%s/%s %s:1\n" %
                                (pdbModelCounter,
                                 str(self.selectedChain1),
                                 alignmentFile1.split("/")[-1]))

            if (self.additionalTargetSequence.get() is True and
                self.inputSequence2.get() is not None and
                self.inputStructureChain.get() is not None):
                f.write("select clear\n")
                f.write("select #%s/%s,%s\n"
                        % (pdbModelCounter,
                           str(self.selectedChain1),
                           str(self.selectedChain2)))

                if self._getOutFastaSequencesFile is not None:
                    alignmentFile2 = self._getOutFastaSequencesFile(self.OUTFILE2)
                    f.write("open %s\n" % alignmentFile2)
                    f.write("sequence disassociate #%s %s\n" %
                            (pdbModelCounter,
                            alignmentFile2.split("/")[-1]))
                    if str(self.selectedModel) != '0':
                        f.write("sequence associate #%s.%s/%s %s:1\n" %
                                (pdbModelCounter,
                                 str(self.selectedModel + 1),
                                 str(self.selectedChain2),
                                 alignmentFile2.split("/")[-1]))
                    else:
                        f.write("sequence associate #%s/%s %s:1\n" %
                                (pdbModelCounter,
                                 str(self.selectedChain2),
                                 alignmentFile2.split("/")[-1]))

        # run the text:
        _chimeraScriptFileName = os.path.abspath(
            self._getTmpPath(chimeraScriptFileName))
        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
            args = " --nogui " + _chimeraScriptFileName
        else:
            args = " " + _chimeraScriptFileName

        f.close()

        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)

        # run in the background
        cwd = os.path.abspath(self._getExtraPath())
        Plugin.runChimeraProgram(Plugin.getProgram(), args, cwd=cwd, extraEnv=getEnvDictionary(self))

    def _validate(self):
        # Check that CLUSTALO or MUSCLE program exists
        errors = super(ChimeraModelFromTemplate, self)._validate()
        if not (self.is_tool(CLUSTALO) or self.is_tool(MUSCLE)):
            errors.append("Clustal-omega and MUSCLE programs missing.\n "
                          "You need at least one of them to run this program.\n"
                          "Please install Clustal-omega and/or MUSCLE:\n"
                          "     sudo apt-get install clustalo\n"
                          "     sudo apt-get install muscle")
        return errors
