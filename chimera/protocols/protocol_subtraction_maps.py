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
from pwem.convert import Ccp4Header
from pwem.objects import Volume
from pwem.objects import Transform
try:
    from pwem.objects import AtomStruct
except ImportError:
    from pwem.objects import PdbFile as AtomStruct

from pyworkflow import VERSION_3_0

from pwem.convert.atom_struct import AtomicStructHandler

from pyworkflow.protocol.params import (PointerParam,
                                        IntParam,
                                        EnumParam,
                                        FloatParam,
                                        BooleanParam,
                                        StringParam, MultiPointerParam)

from pwem.constants import (SYM_DIHEDRAL_X)
from ..constants import (CHIMERA_I222, CHIMERA_SYM_NAME)
from ..convert import CHIMERA_LIST
from pwem.protocols import EMProtocol
from .protocol_base import createScriptFile
from pwem.viewers.viewer_chimera import (Chimera,
                                         sessionFile,
                                         chimeraMapTemplateFileName,
                                         chimeraScriptFileName,
                                         chimeraPdbTemplateFileName)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .. import Plugin

class ChimeraSubtractionMaps(EMProtocol):
    """Protocol to subtract two volumes.
        One of these volumes can be derived from an atomic structure.
        Execute command *scipionwrite [model #n]* from command line in order
        to transfer the selected
        pdb to scipion. Default value is model=#0,
        model refers to the pdb file."""
    _label = 'map subtraction'
    _program = ""
    _version = VERSION_3_0

    MAP_OPTIONS = ['3D map', 'atomic structure']
    CHIMERA_FILTERS = ['Gaussian', 'Fourier Transform']

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=False):

        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input 3D Map',
                      important=True,
                      help="difference 3D map = minuend − subtrahend"
                           "input here the minuend")
        form.addParam('mapOrModel', EnumParam,
                     choices=self.MAP_OPTIONS,
                     display=EnumParam.DISPLAY_HLIST,
                     default=0, label='Subtraction of',
                     help="difference 3D Map = minuend − subtrahend. "
                          "Subtrahend 3D map may be provided by the user (choose '3D map') "
                          "or created from an atomic coordinates file (choose 'atomic structure'"
                          " If 3D Map is chosen, the sampling rate of minuend should be"
                          " equal to sampling rate of subtrahend")
        form.addParam('inputVolume2', PointerParam, pointerClass="Volume",
                         condition=('mapOrModel==%d ' % 0),
                         important=True, allowsNull=True,
                         label='Map to subtract (subtrahend)',
                         help="Map that has to be subtracted from the minuend 3D map")
        form.addParam('resolution', FloatParam,
                         condition=('mapOrModel==%d' % 1),
                         label='Map resolution (A):',
                         help="Thershold that will be used to low pass the 3D map created "
                              "from the atomic structure. We recommend half the value of the"
                              " resolution obtained by FSC")
        form.addParam('pdbFileToBeRefined', PointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      condition=('mapOrModel==%d ' % 1),
                      important=True,
                      label='Atomic structure',
                      help="Atomic structure to derive "
                       "a 3D map that will be subtracted from the minuend map.")
        form.addParam('selectChain', BooleanParam,
                         condition=('mapOrModel==%d' % 1),
                         label="Select a specific chain?",
                         default=False,
                         help="Select 'Yes' if you want to generate the map"
                              " from a specific chain of the atomic structure.\n")
        form.addParam('selectStructureChain', StringParam,
                         condition=('mapOrModel==%d and selectChain==True' % 1),
                         label="Chain of the atomic structure", important=True,
                         help="Select a particular chain of the atomic "
                              "structure.")
        form.addParam('removeResidues', BooleanParam,
                         condition=('mapOrModel==%d' % 1),
                         label="Remove residues from the atomic structure?",
                         default=False,
                         help="Select 'Yes' to remove a certain "
                              "number of residues of one of the chains "
                              "of the atomic structure. these removed residues "
                              "might help you to establish a control of "
                              "appropriate levels of map density.\n")
        form.addParam('inputStructureChain', StringParam,
                         condition=(('mapOrModel==%d and '
                                    'removeResidues==True and '
                                    'selectChain==False') % 1),
                         label="Chain ", important=True,
                         help="Select a particular chain of the atomic "
                              "structure.")
        form.addParam('firstResidueToRemove', StringParam,
                         condition=(('mapOrModel==%d and '
                                    'removeResidues==True') % 1),
                         label="First residue to remove ", default=1,
                         important=True,
                         help="Select the first residue of the selected atomic "
                              "structure chain that you want to remove.")
        form.addParam('lastResidueToRemove', StringParam,
                         condition=(('mapOrModel==%d and '
                                    'removeResidues==True') % 1),
                         label="Last residue to remove ", default=2,
                         important=True,
                         help="Select the last residue of the selected atomic "
                              "structure chain that you want to remove.")
        form.addParam('applySymmetry', BooleanParam,
                      condition=('mapOrModel==%d' % 1),
                      label="Apply symmetry to the atomic structure:",
                      default=False,
                      help="'Symmetry = Yes' indicates that symmetry will be applied. "
                           "This option is recommended is the atomic structure "
                           "corresponds to the asymmetrical unit and you want to "
                           "regenerate the structure of the whole map.\n'Symmetry = "
                           "No' indicates that symmetry will not be applied, and then  "
                           "the map derived from the atomic structure involves the input"
                           "atomic structure only.\n")
        form.addParam('symmetryGroup', EnumParam,
                      condition=(('mapOrModel==%d and applySymmetry==True') % 1),
                      choices=CHIMERA_LIST,
                      default=CHIMERA_I222,
                      important=True,
                      label="Symmetry",
                      help="https://scipion-em.github.io/docs/release-2.0.0/docs/developer/symmetries.html?highlight=symmetry"
                           "Symmetry for a description of the symmetry groups "
                           "format in CHIMERA.\n"
                           "If no symmetry is present, use _c1_."
                           'More information: \n'
                           'https://www.cgl.ucsf.edu/chimera/current/docs/UsersGuide/midas/sym.html'
                      )
        form.addParam('symmetryOrder', IntParam, default=1,
                      condition=(('mapOrModel==%d' % 1) and
                                 'applySymmetry==True' and
                                 ('symmetryGroup<=%d' % SYM_DIHEDRAL_X)),
                      label='Symmetry Order',
                      help='Select the order of cyclic or dihedral symmetry.')
        form.addParam('rangeDist', IntParam, default=100,
                      condition=(('mapOrModel==%d and applySymmetry==True') % 1),
                      label='Range of distance',
                      help="This value allows to generate copies with centers within a certain range "
                           "of distance of the center of the original molecule model. A models's center "
                           "is defined as the center of its bounding box.")
        form.addParam('inputPdbFiles', MultiPointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      label='Other atomic structures',
                      help="In case you need to load more PDBx/mmCIF files, "
                           "you can load them here. This files will NOT be used "
                           "to create the subtrahend")
        form.addParam('selectAreaMap', BooleanParam,
                      condition=('mapOrModel==%d' % 1),
                      label="Map fraction around the atomic structure?",
                      default=False,
                      help="Select 'Yes' if you want to limit the map"
                           " to a certain radius around the atomic structure.\n")
        form.addParam('radius', IntParam, default=15,
                      condition=(('mapOrModel==%d' % 1) and
                                 'selectAreaMap==True'),
                      label="Atom radius (Angstroms)",
                      help="Set the radius to select values of grid points "
                           "farther than that radius from any atom.")
        form.addParam("filterToApplyToDiffMap", EnumParam,
                      expertLevel=LEVEL_ADVANCED,
                      choices=self.CHIMERA_FILTERS,
                      display=EnumParam.DISPLAY_HLIST,
                      label="Filter to apply to the differential map",
                      default=0,
                      help="Choose the filter to clean the background noise "
                           "of the differential map.")
        form.addParam("widthFilter", FloatParam,
                      condition='filterToApplyToDiffMap',
                      expertLevel=LEVEL_ADVANCED,
                      label="Gaussian filter width",
                      default=1.5,
                      help="Set the width of the Gaussian filter.")
        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="Add extra commands in cmd file. Use for testing")
        form.addSection(label='Help')
        form.addLine("Step 1:\nIn the sequence window your target "
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

    def _insertAllSteps(self):
        self._insertFunctionStep('prerequisitesStep')
        self._insertFunctionStep('runChimeraStep')
        self._insertFunctionStep('createOutput')


    def prerequisitesStep(self):
        # get the reference map
        self.vol = self.inputVolume.get()
        self.fnVolName = os.path.abspath(self.vol.getFileName())
        print("Reference map:\n %s\n" % self.fnVolName)

        # clean form
        if self.mapOrModel.get() == 0: # 0 -> map
            self.subVol = self.inputVolume2.get()
            self.subVolName = os.path.abspath(self.subVol.getFileName())
            print("Map to subtract:\n %s\n" % self.subVolName)
        else: # 1 -> ATOMIC STRUCT
            self.atomStruct = self.pdbFileToBeRefined.get()
            self.atomStructName = os.path.abspath(self.atomStruct.getFileName())
            print("Map to subtract generated from the atomic structure:\n %s\n"
                  % self.atomStructName)

    def runChimeraStep(self):
        # building script file including the coordinate axes and the input
        # volume with samplingRate and Origin information
        f = open(self._getTmpPath(chimeraScriptFileName), "w")
        f.write("from chimera import runCommand\n")

        # create coherent header
        createScriptFile(1,  # model id pdb
                         1,  # model id 3D map
                         self._getExtraPath(chimeraPdbTemplateFileName),
                         self._getExtraPath(chimeraMapTemplateFileName),
                         f,
                         self._getExtraPath(sessionFile),
                         )

        # building coordinate axes
        dim = self.vol.getDim()[0]
        sampling = self.vol.getSamplingRate()
        bildFileName = os.path.abspath(self._getTmpPath("axis_input.bild"))
        Chimera.createCoordinateAxisFile(dim, bildFileName=bildFileName,
                                         sampling=sampling)
        # input volume
        chimeraModelId = 0
        f.write("runCommand('open %s')\n" % (bildFileName))
        f.write("runCommand('cofr 0,0,0')\n")  # set center of coordinates
        # origin coordinates
        chimeraModelMapM = chimeraModelId + 1 # 1, Minuend
        f.write("runCommand('open %s')\n" % self.fnVolName)
        f.write("runCommand('volume #%d style surface voxelSize %f')\n"
                % (chimeraModelMapM, sampling))
        x, y, z = self.vol.getShiftsFromOrigin()
        f.write("runCommand('volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                % (chimeraModelMapM, x, y, z))
        # input vol with its origin coordinates
        chimeraModelMapS = chimeraModelMapM + 1  # 2 Subtrahend
        if self.mapOrModel == 0:
            f.write("runCommand('open %s')\n" %
                    (self.subVolName))
            f.write("runCommand('volume #%d style surface voxelSize %f')\n"
                    % (chimeraModelMapS, sampling))
            x, y, z = self.subVol.getShiftsFromOrigin()
            f.write("runCommand('volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                        % (chimeraModelMapS, x, y, z))
            chimeraModelMapDiff = chimeraModelMapS + 1  # 3 Diff map
            f.write("runCommand('vop subtract #%d #%d modelId #%d "
                    "minRMS true onGrid #%d')\n"
                    % (chimeraModelMapM, chimeraModelMapS,
                       chimeraModelMapDiff, chimeraModelId + 1))
            f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                    "prefix difference_')\n"
                    % (chimeraModelMapDiff, chimeraModelMapM))
        else:
            f.write("runCommand('open %s')\n"  % (self.atomStructName))
            chimeraModelAtomStruct = chimeraModelMapM + 1
            if self.selectChain == True:
                # model and chain selected
                if self.selectStructureChain.get() is not None:
                    chain = self.selectStructureChain.get()
                    self.selectedModel = chain.split(',')[0].split(':')[1].strip()
                    chimeraModelAtomStruct = int(chimeraModelAtomStruct +
                                                 int(self.selectedModel))
                    self.selectedChain = \
                        chain.split(',')[1].split(':')[1].strip().split('"')[1]
                    print("Selected chain: %s from model: %s from structure: %s" \
                        % (self.selectedChain, self.selectedModel,
                            os.path.basename(self.atomStructName)))
                    f.write("runCommand('sel #%d:.%s')\n"
                            % (chimeraModelAtomStruct, self.selectedChain))
                    f.write("runCommand('write format pdb selected relative %d #%d "
                            "/tmp/chain.pdb')\n"
                            % (chimeraModelId, chimeraModelAtomStruct))
                    f.write("runCommand('open /tmp/chain.pdb')\n")
                    chimeraModelAtomStructChain = chimeraModelAtomStruct + 1
                    f.write("runCommand('scipionwrite model #%d refmodel #%d "
                            "prefix chain_%s_')\n"
                            % (chimeraModelAtomStructChain, chimeraModelMapM,
                               self.selectedChain))
                    if self.applySymmetry == True:
                        if self.symmetryGroup.get() is not None:
                            sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                            if sym == "Cn" and self.symmetryOrder != 1:
                                print("self.symmetryOrder: ", self.symmetryOrder)
                                print("self.rangeDist: ", self.rangeDist)
                                f.write("runCommand('sym #%d group C%d range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           self.symmetryOrder, self.rangeDist))
                            elif sym == "Dn" and self.symOrder != 1:
                                f.write("runCommand('sym #%d group d%d range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           self.symmetryOrder, self.rangeDist))
                            elif sym == "T222" or self.sym == "TZ3":
                                f.write("runCommand('sym #%d group t,%s range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           sym[1:], self.rangeDist))
                            elif sym == "O":
                                f.write("runCommand('sym #%d group O range %d')\n"
                                        % chimeraModelAtomStructChain, self.rangeDist)
                            elif sym == "I222" or sym == "I222r" or sym == "In25" or \
                                    sym == "In25r" or sym == "I2n3" or sym == "I2n3r" or \
                                    sym == "I2n5" or sym == "I2n5r":
                                f.write("runCommand('sym #%d group i,%s range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           sym[1:], self.rangeDist))
                            chimeraModelAtomStructChainSym = chimeraModelAtomStructChain + 2
                            f.write("runCommand('combine #%d- modelId #%d')\n"
                                    % (chimeraModelAtomStructChain, chimeraModelAtomStructChainSym))
                            if self.selectAreaMap == True:
                                chimeraModelIdZone = chimeraModelAtomStructChainSym + 1
                                f.write("runCommand('vop zone #%d #%d %d modelId #%d')\n"
                                        % (chimeraModelMapM, chimeraModelAtomStructChainSym,
                                           self.radius, chimeraModelIdZone))
                                f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                                        "prefix zone_')\n" % (chimeraModelIdZone,
                                                              chimeraModelMapM))
                                f.write("runCommand('close #%d')\n" % (chimeraModelAtomStructChain + 1))
                                f.write("runCommand('close #%d')\n" % (chimeraModelAtomStructChainSym))
                    else:
                        if self.selectAreaMap == True:
                            chimeraModelIdZone = chimeraModelAtomStructChain + 1
                            f.write("runCommand('vop zone #%d #%d %d modelId #%d')\n"
                                    % (chimeraModelMapM, chimeraModelAtomStructChain,
                                    self.radius, chimeraModelIdZone))

                    chimeraModelAtomStructChainMap = chimeraModelIdZone + 1
                    chimeraModelMapDiff = chimeraModelAtomStructChainMap + 1

                    if self.removeResidues == True:
                        if (self.firstResidueToRemove.get() is not None and
                                self.lastResidueToRemove.get() is not None):
                            firstResidue = self.firstResidueToRemove.get().\
                            split(":")[1].split(",")[0].strip()
                            lastResidue = self.lastResidueToRemove.get(). \
                                split(":")[1].split(",")[0].strip()
                            print("firstResidue: ", firstResidue)
                            print("lastResidue: ", lastResidue)
                            f.write("runCommand('select #%d:%d-%d.%s')\n"
                                    % (chimeraModelAtomStructChain,
                                       int(firstResidue), int(lastResidue),
                                       self.selectedChain))
                            f.write("runCommand('del sel')\n")
                    if self.applySymmetry == True:
                        if self.symmetryGroup.get() is not None:
                            sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                            if sym == "Cn" and self.symmetryOrder != 1:
                                f.write("runCommand('sym #%d group C%d range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           self.symmetryOrder, self.rangeDist))
                            elif sym == "Dn" and self.symOrder != 1:
                                f.write("runCommand('sym #%d group d%d range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           self.symmetryOrder, self.rangeDist))
                            elif sym == "T222" or self.sym == "TZ3":
                                f.write("runCommand('sym #%d group t,%s range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           sym[1:], self.rangeDist))
                            elif sym == "O":
                                f.write("runCommand('sym #%d group O range %d')\n"
                                        % chimeraModelAtomStructChain, self.rangeDist)
                            elif sym == "I222" or sym == "I222r" or sym == "In25" or \
                                    sym == "In25r" or sym == "I2n3" or sym == "I2n3r" or \
                                    sym == "I2n5" or sym == "I2n5r":
                                f.write("runCommand('sym #%d group i,%s range %d')\n"
                                        % (chimeraModelAtomStructChain,
                                           sym[1:], self.rangeDist))
                            chimeraModelAtomStructChainSym = chimeraModelAtomStructChain + 2
                            f.write("runCommand('combine #%d- modelId #%d')\n"
                                    % (chimeraModelAtomStructChain, chimeraModelAtomStructChainSym))

                            f.write("runCommand("
                                    "'molmap #%d %0.3f gridSpacing %0.2f modelId #%d')\n"
                                    % (chimeraModelAtomStructChainSym, self.resolution, sampling,
                                       chimeraModelAtomStructChainMap))
                    else:
                        f.write("runCommand("
                                "'molmap #%d %0.3f gridSpacing %0.2f modelId #%d')\n"
                                % (chimeraModelAtomStructChain, self.resolution, sampling,
                                   chimeraModelAtomStructChainMap))
                    f.write("runCommand('scipionwrite model #%d refmodel #%d "
                            "prefix molmap_chain%s_')\n"
                            % (chimeraModelAtomStructChainMap, chimeraModelMapM,
                               self.selectedChain))

                    if self.selectAreaMap == True:
                        f.write("runCommand('vop subtract #%d #%d modelId #%d "
                            "minRMS true onGrid #%d')\n"
                            % (chimeraModelIdZone, chimeraModelAtomStructChainMap,
                                chimeraModelMapDiff, chimeraModelMapM))
                    else:
                        f.write("runCommand('vop subtract #%d #%d modelId #%d "
                                "minRMS true onGrid #%d')\n"
                                % (chimeraModelMapM, chimeraModelAtomStructChainMap,
                                   chimeraModelMapDiff, chimeraModelMapM))
                    f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                            "prefix difference_')\n"
                            % (chimeraModelMapDiff, chimeraModelMapM))
        #             else:
        #                 chimeraModelMapDiff = chimeraModelAtomStructChainMap + 1  # 3 Diff map
        #                 f.write("runCommand('vop subtract #%d #%d modelId #%d "
        #                         "minRMS true onGrid #%d')\n"
        #                         % (chimeraModelMapM, chimeraModelAtomStructChainMap,
        #                            chimeraModelMapDiff, chimeraModelMapM))
        #                 f.write("runCommand('scipionwrite model #%d refmodel #%d " \
        #                         "prefix difference_')\n"
        #                         % (chimeraModelMapDiff, chimeraModelMapM))
        #
        #
        #
        #         # if self.removeResidues == True:
        #         #     pass
        #         # else:
        #         #     f.write("runCommand('open %s')\n" )
        #         #
        #         #
        #         # # f.write("runCommand('')")
        #     else:
        #         if self.selectAreaMap == True:
        #             chimeraModelIdZone = chimeraModelAtomStruct + 1
        #             f.write("runCommand('vop zone #%d #%d %d modelId #%d')\n"
        #                     % (chimeraModelMapM, chimeraModelAtomStruct,
        #                        self.radius, chimeraModelIdZone))
        #             f.write("runCommand('scipionwrite model #%d refmodel #%d " \
        #                     "prefix zone_')\n"
        #                     % (chimeraModelIdZone, chimeraModelMapM))
        #
        #     if self.removeResidues == True:
        #         pass
        #
        #

        chimeraModelMapDiffFil = chimeraModelMapDiff + 1
        if self.filterToApplyToDiffMap.get() == 0:
            f.write("runCommand('vop gaussian #%d sd %0.3f')\n"
                    % (chimeraModelMapDiff, self.widthFilter.get()))
        else:
            f.write("rc('vop laplacian #%d')\n"
                    % (chimeraModelMapDiff))
        f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                "prefix filtered_')\n"
                % (chimeraModelMapDiffFil, chimeraModelMapM))

        # run the text:
        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
            args = " --nogui --script " + self._getTmpPath(
                            chimeraScriptFileName)
        else:
            args = " --script " + self._getTmpPath(chimeraScriptFileName)
        f.close()

        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)

        # run in the background
        Chimera.runProgram(Plugin.getProgram(), args)

    def createOutput(self):

        """ Register outputs.
        """

        # Check vol and pdb files
        directory = self._getExtraPath()
        counter = 1
        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".mrc"):
                volFileName = os.path.join(directory, filename)
                vol = Volume()
                vol.setFileName(volFileName)

                # fix mrc header
                ccp4header = Ccp4Header(volFileName, readHeader=True)
                sampling = ccp4header.computeSampling()
                origin = Transform()
                shifts = ccp4header.getOrigin()
                origin.setShiftsTuple(shifts)
                vol.setOrigin(origin)
                vol.setSamplingRate(sampling)
                keyword = filename.split(".mrc")[0]
                # counter += 1
                kwargs = {keyword: vol}
                self._defineOutputs(**kwargs)

            if filename.endswith(".pdb") or filename.endswith(".cif"):
                path = os.path.join(directory, filename)
                pdb = AtomStruct()
                pdb.setFileName(path)
                keyword = filename.split(".pdb")[0].replace(".","_")
                #counter += 1
                kwargs = {keyword: pdb}
                self._defineOutputs(**kwargs)


    def _validate(self):
        errors = super(ChimeraSubtractionMaps, self)._validate()

        # Check that the input volume exists
        if self.inputVolume.get() is None:
            errors.append("Error: You should provide a map.\n")

        # Check that an additional map or an atomic structure exist
        if (not self.inputVolume2.get() and
                not self.pdbFileToBeRefined.get()):
            errors.append("Error: You should provide an additional "
                          "map or an atomic structure.\n")

        return errors

