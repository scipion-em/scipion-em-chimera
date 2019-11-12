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
from pyworkflow.em import *
from pyworkflow import VERSION_1_2
from chimera.protocols import ChimeraProtBase

# Horrible hack to release this plugin before scipion next version.
# TODO: remove when possible
from pyworkflow import LAST_VERSION, VERSION_2_0
if LAST_VERSION == VERSION_2_0 :
    from chimera.atom_struct import AtomicStructHandler
else:
    from pyworkflow.em.convert.atom_struct import AtomicStructHandler

from pyworkflow.protocol.params import (PointerParam,
                                        StringParam,
                                        MultiPointerParam)
from pyworkflow.em.convert.sequence import (SequenceHandler,
                                            saveFileSequencesToAlign,
                                            alignClustalSequences,
                                            alignBioPairwise2Sequences,
                                            alignMuscleSequences)
from collections import OrderedDict
from chimera.constants import CLUSTALO, MUSCLE

class ChimeraModelFromTemplate(ChimeraProtBase):
    """Protocol to model three-dimensional structures of proteins using Modeller.
        Execute command *scipionwrite [model #n]* from command line in order
        to transfer the selected
        pdb to scipion. Default value is model=#0,
        model refers to the pdb file."""
    _label = 'model from template'
    _program = ""
    _version = VERSION_1_2

    SEQUENCEFILENAME = '_sequence.fasta'
    INFILE = "unaligned.fasta"
    OUTFILE = "aligned.fasta"
    TWOSEQUENCES = 0
    MULTIPLESEQUENCES = 1
    ProgramToAlign1 = ['Bio.pairwise2', 'Clustal Omega', 'MUSCLE']
    ProgramToAlign2 = ['Clustal Omega', 'MUSCLE']

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
        # hide inputPdbFiles
        param = form.getParam('inputPdbFiles')
        param.condition.set('False')
        param.allowsNull.set('True')
        section = formBase.getSection('Input')
        section.addParam('inputStructureChain', StringParam,
                       label="Chain ", allowsNull=True, important=True,
                       help="Select a particular chain of the atomic "
                            "structure.")
        section.addParam('inputSequence', PointerParam, pointerClass="Sequence",
                       label='Target sequence', allowsNull=True,
                       important=True,
                       help="Input the aminoacid sequence to align with the "
                            "structure template sequence.")
        section.addParam('additionalSequencesToAlign', params.BooleanParam,
                         default=False, label='Additional sequences to align?',
                         help='Select YES if you want to add some more '
                              'sequences to accomplish the alignment. This '
                              'option is recommendable when the first two '
                              'sequences are not very similar.')
        section.addParam('inputSequencesToAlign', MultiPointerParam,
                         pointerClass="Sequence", allowsNull=True,
                         condition='additionalSequencesToAlign == True',
                         label='Other sequences to align',
                         help="In case you need to load more sequences to "
                              "align, you can load them here.")
        section.addParam('inputProgramToAlign1', params.EnumParam,
                         choices=self.ProgramToAlign1,
                         label="Alignment tool for two sequences:", default=0,
                         condition='additionalSequencesToAlign == False',
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
        section.addParam('inputProgramToAlign2', params.EnumParam,
                         choices=self.ProgramToAlign2,
                         label="Multiple alignment tool:", default=0,
                         condition='additionalSequencesToAlign == True',
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
        and write the Modeller license key supplied (Academic users can 
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
        # get pdb sequence
        import json
        chainIdDict = json.loads(self.inputStructureChain.get())

        # read PDB
        self.structureHandler = AtomicStructHandler()
        fileName = os.path.abspath(self.pdbFileToBeRefined.get(
            ).getFileName())
        self.structureHandler.read(fileName)

        # get sequence of structure chain with id chainId (selected by the user)
        self.selectedModel = chainIdDict['model']
        self.selectedChain = chainIdDict['chain']
        # self.selectedModel = chainId.split(',')[0].split(':')[1].strip()
        # self.selectedChain = chainId.split(',')[1].split(':')[1].strip()
        print "Selected chain: %s from model: %s from structure: %s" \
              % (self.selectedChain, self.selectedModel,
                 os.path.basename(fileName))

        # Bio.Seq.Seq object
        structureSeq = self.structureHandler.getSequenceFromChain(
            self.selectedModel, self.selectedChain)

        # obtain a seqID for our PDB sequence
        structSeqID = self.structureHandler.getFullID(self.selectedModel,
                                                self.selectedChain)
        # END PDB sequence

        # start user imported target sequence
        # get target sequence imported by the user
        userSeq = self.inputSequence.get()  # SEQ object from Scipion
        targetSeqID = userSeq.getId()  # ID associated to SEQ object (str)
        userSequence = userSeq.getSequence()  # sequence associated to
                                                   # that SEQ object (str)
        # transformation of this sequence (str) in a Bio.Seq.Seq object:
        seqHandler = SequenceHandler(userSequence,
                                     isAminoacid=userSeq.getIsAminoacids())
        targetSeq = seqHandler._sequence # Bio.Seq.Seq object

        # creation of Dic of IDs and sequences
        SeqDic = OrderedDict()
        SeqDic[structSeqID] = structureSeq
        SeqDic[targetSeqID] = targetSeq

        # align sequences and save them to disk, -this will be chimera input-
        # get all sequences in a fasta file
        inFile = self._getInFastaSequencesFile()
        outFile = self._getOutFastaSequencesFile()

        # get the alignment of sequences
        if self.additionalSequencesToAlign.get() == False:
            saveFileSequencesToAlign(SeqDic, inFile)
            self.inputSequencesToAlign = None
            if self.inputProgramToAlign1.get() == \
                    self.ProgramToAlign1.index('Bio.pairwise2'):
            # Only the two first sequences will be included in the alignment
                self.alignment = alignBioPairwise2Sequences(
                    structSeqID, structureSeq,
                    targetSeqID, targetSeq,
                    outFile)
            else:
                # All the sequences will be included in the alignment
                if self.inputProgramToAlign1.get() == \
                        self.ProgramToAlign1.index('Clustal Omega'):
                    cline = alignClustalSequences(inFile, outFile)
                else:
                    cline = alignMuscleSequences(inFile, outFile)
                args = ''
                self.runJob(cline, args)
        else:
            # if there are additional sequences imported by the user
            if self.inputSequencesToAlign is not None:
                for seq in self.inputSequencesToAlign:
                    seq = seq.get()
                    ID = seq.getId()
                    sequence = seq.getSequence()
                    seqHandler = SequenceHandler(sequence,
                                                 isAminoacid=seq.getIsAminoacids())
                    otherSeq = seqHandler._sequence  # Bio.Seq.Seq object
                    SeqDic[ID] = otherSeq

            # align sequences and save them to disk, -this will be chimera input-
            # get all sequences in a fasta file
            #inFile = self._getInFastaSequencesFile()
            saveFileSequencesToAlign(SeqDic, inFile)
            #outFile = self._getOutFastaSequencesFile()

            # All the sequences will be included in the alignment
            if self.inputProgramToAlign2.get() == self.ProgramToAlign2.index(
                    'Clustal Omega'):
                cline = alignClustalSequences(inFile, outFile)
            else:
                cline = alignMuscleSequences(inFile, outFile)
            args = ''
            self.runJob(cline, args)


    def _getInFastaSequencesFile(self):
        INFILENAME = self._getTmpPath(self.INFILE)
        return os.path.abspath(INFILENAME)

    def _getOutFastaSequencesFile(self):
        OUTFILENAME = self._getExtraPath(self.OUTFILE)
        return os.path.abspath(OUTFILENAME)

    #def validate(self):
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
