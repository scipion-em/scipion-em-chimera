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
from pyworkflow import VERSION_1_2
from pyworkflow.utils import copyFile

from . import ChimeraProtBase

from pwem.convert.atom_struct import AtomicStructHandler

from pyworkflow.protocol.params import (PointerParam,
                                        StringParam,
                                        MultiPointerParam,
                                        BooleanParam,
                                        EnumParam,
                                        PathParam)
from pwem.convert.sequence import (SequenceHandler,
                                   saveFileSequencesToAlign,
                                   alignClustalSequences,
                                   alignBioPairwise2Sequences,
                                   alignMuscleSequences)
from collections import OrderedDict
from ..constants import CLUSTALO, MUSCLE


class ChimeraModelFromTemplate(ChimeraProtBase):
    """Protocol to model three-dimensional structures of proteins using Modeller.
        Execute command *scipionwrite [model #n]* from command line in order
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

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=False):
        formBase = super(ChimeraModelFromTemplate, self)._defineParams(form,
                                                                       doHelp=True)
        param = form.getParam('pdbFileToBeRefined')
        param.label.set('Atomic structure used as template')
        param.help.set("PDBx/mmCIF file template used as basic atomic "
                       "structure to model your specific sequence.")
        param = form.getParam('inputVolume')
        param.condition.set('False')
        param = form.getParam('inputVolumes')
        param.condition.set('False')
        # hide inputPdbFiles
        param = form.getParam('inputPdbFiles')
        param.condition.set('False')
        param.allowsNull.set('True')
        section = formBase.getSection('Input')
        section.addParam('inputStructureChain', StringParam,
                         label="Chain ", allowsNull=True, important=True,
                         help="Select a particular chain of the atomic "
                              "structure.")
        section.addParam('inputSequence1', PointerParam, pointerClass="Sequence",
                         label='Target sequence', allowsNull=True,
                         important=True,
                         help="Input the aminoacid sequence to align with the "
                              "structure template sequence.")
        section.addParam('optionForAligning1', EnumParam,
                         choices=self.OptionForAligning,
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
                         condition='optionForAligning1 == 2',
                         label='Sequence alignment input',
                         help="Input your own sequence alignment.\n"
                              "ChimeraX allowed formats accessible here: "
                              "https://www.cgl.ucsf.edu/chimerax/docs/user/commands/open.html#sequence ")
        section.addParam('inputSequencesToAlign1', MultiPointerParam,
                         pointerClass="Sequence", allowsNull=True,
                         condition='optionForAligning1 == 1',
                         label='Other sequences to align',
                         help="In case you need to load more sequences to "
                              "align, you can load them here.")
        section.addParam('inputProgramToAlign1_1', EnumParam,
                         choices=self.ProgramToAlign1,
                         label="Alignment tool for two sequences:", default=0,
                         condition='optionForAligning1 == 0',
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
                         condition='optionForAligning1 == 1',
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
                         default=False, label='Additional target sequence to include?',
                         help='Select YES if you want to add an additional '
                              'target sequence to model according a different '
                              'chain of the structure template. This '
                              'option is recommendable when you want to model '
                              'the two interacting elements of a particular complex'
                              ' at the same time.')
        section.addParam('selectStructureChain', StringParam,
                         condition='additionalTargetSequence == True',
                         label="Chain ", allowsNull=True, important=True,
                         help="Select a particular chain of the atomic "
                              "structure.")
        section.addParam('inputSequence2', PointerParam, pointerClass="Sequence",
                         condition='additionalTargetSequence == True',
                         label='Target sequence', allowsNull=True,
                         important=True,
                         help="Input the aminoacid sequence to align with the "
                              "structure template sequence.")
        section.addParam('optionForAligning2', EnumParam,
                         choices=self.OptionForAligning,
                         condition='additionalTargetSequence == True',
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
                         condition='optionForAligning2 == 2 and '
                                   'additionalTargetSequence == True',
                         label='Sequence alignment input',
                         help="Input your own sequence alignment.\n"
                              "ChimeraX allowed formats accessible here: "
                              "https://www.cgl.ucsf.edu/chimerax/docs/user/commands/open.html#sequence ")
        section.addParam('inputSequencesToAlign2', MultiPointerParam,
                         pointerClass="Sequence", allowsNull=True,
                         condition='optionForAligning2 == 1 and '
                                   'additionalTargetSequence == True',
                         label='Other sequences to align',
                         help="In case you need to load more sequences to "
                              "align, you can load them here.")
        section.addParam('inputProgramToAlign1_2', EnumParam,
                         choices=self.ProgramToAlign1,
                         label="Alignment tool for two sequences:", default=0,
                         condition='optionForAligning2 == 0 and '
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
                         condition='optionForAligning2 == 1 and '
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
        formBase.addLine("Step 1:\nIn the sequence window your target "
                         "sequence (and other additional sequences that you "
                         "want to use in  the alignment) will appear aligned to "
                         "the template's sequence. Select in the sequence window "
                         "menu:\nStructure -> Modeller (homology)...;\nA new  "
                         "window for Comparative Modeling with Modeller will "
                         "appear. Select your specific sequence as the sequence "
                         "to be modeled (target), and the input atomic structure"
                         + '''
        used as template for modeling. Select Run Modeller via web service 
        and write the Modeller license key supplied (Academic user can 
        register free of charge to receive a license key). Finally, press OK.
        \nWAITING TIME: (you may see the status of your job in chimera main 
        window, lower left corner.)\n\nStep 2:\nWhen the process finished, 
        5 models will 
        be automatically superimposed onto the template and model scores
        will appear in Modeller Results window. In Chimera main menu -> 
        Favorites -> Model panel will show you: #0 (coordinate axes); #1 (
        template); #2.1 to 2.5 (models).Choose the one you like the best, 
        for example model #2.1. To save it in Scipion, we need to change the 
        model ID. In Chimera main menu: Favorites -> Command Line, write 
        *combine #2.1 model #3 close t*. Then, you will see in Model panel 
        that selected model #2.1 renamed to combination with ID #3. Save it 
        as first guess in Scipion by executing the Chimera command 
        *scipionwrite [model #n]*. In our example *scipionwrite model #3*.\n 
        When you use the command line scipionwrite, the Chimera session will 
        be saved by default. Additionally, you can save the Chimera session 
        whenever you want by executing the command *scipionss*. You will be 
        able to restore the saved session by using the protocol chimera 
        restore session (SCIPION menu: Tools/Calculators/chimera restore 
        session). Once you have save your favorite model you can press 
        Quit in the Modeller Results window.''')

    # --------------------------- INSERT steps functions --------------------

    def prerequisitesStep(self):

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

    def prePreRequisites(self, fileName, chainIdDict, userSeq, inFile, \
                         outFile, addSeq, yourAlignment, inputSeqAlign, \
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

    # def validate(self):
    #    super(ChimeraModelFromTemplate, self).validate()
    #    # TODO check if clustal/muscle exists
    #    #TODO; change help ro installation pages instead of apt-get

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
