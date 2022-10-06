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

import json

from pwem import *
from pwem.convert import Ccp4Header
from pwem.objects import Volume
from pwem.objects import Transform
#try:
from pwem.objects import AtomStruct
#except ImportError:
#    from pwem.objects import PdbFile as AtomStruct

from pyworkflow import VERSION_3_0

from pyworkflow.protocol.params import (PointerParam,
                                        IntParam,
                                        EnumParam,
                                        FloatParam,
                                        BooleanParam,
                                        StringParam, MultiPointerParam, Float)

from pwem.constants import (SYM_DIHEDRAL_X)
from ..constants import (CHIMERA_SYM_NAME, CHIMERA_I222r)
from ..convert import CHIMERA_LIST
from pwem.protocols import EMProtocol
from pwem.viewers.viewer_chimera import Chimera, chimeraScriptFileName
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from chimera import Plugin
from pyworkflow.utils.properties import Message

from chimera.utils import getEnvDictionary


class ChimeraSubtractionMaps(EMProtocol):
    """Protocol to subtract two volumes.
        One of these volumes can be derived from an atomic structure.
        Execute command *scipionwrite #n [prefix stringAddedToFilename]*
        from command line in order to transfer the generated maps and models to scipion.
        In addition to maps and models that the protocol saves by default,
        the user can generate and save some others"""
    _label = 'map subtraction'
    _program = ""
    _version = VERSION_3_0

    @classmethod
    def getClassPackageName(cls):
        return "chimerax"

    PROTOCOL_OPTIONS = ['Subtraction', 'Mask']
    MAP_OPTIONS = ['3D map', 'atomic structure']
    CHIMERA_FILTERS = ['Gaussian', 'Fourier Transform']
    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=False):

        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input 3D Map',
                      important=True,
                      help="Difference 3D map = minuend − subtrahend.\n"
                           "Input here the minuend of the subtraction.")
        form.addParam('subtractOrMask', EnumParam,
                      choices=self.PROTOCOL_OPTIONS,
                      display=EnumParam.DISPLAY_HLIST,
                      default=0,
                      label='Select the operation to perform',
                      help='You can select "Subtract" to get the result '
                           'minuend − subtrahend, or "Mask" to mask the '
                           'minuend.\nThis mask is created with all '
                           'those points that belong to the subtrahend '
                           'and are greater than the level (0 level is '
                           'not supplied).')
        form.addParam('level', FloatParam,
                      expertLevel=LEVEL_ADVANCED,
                      condition=('subtractOrMask==%d ' % 1),
                      default=0.001,
                      allowsNull=True,
                      label='Contour level (subtrahend)',
                      help='Difference 3D map = minuend − subtrahend.\n"'
                           'Calculation are made for those voxels '
                           'inside a region created by this contour level\n'
                           'Empty value -> chimera computes the level.')
        form.addParam('mapOrModel', EnumParam,
                     choices=self.MAP_OPTIONS,
                     display=EnumParam.DISPLAY_HLIST,
                     default=0, label='Subtraction/Mask of',
                     help="Difference 3D map = minuend − subtrahend.\n"
                          "Subtrahend 3D map may be provided by the user "
                          "(choose '3D map') "
                          "or created from an atomic coordinates file "
                          "(choose 'atomic structure').\n"
                          " If 3D Map is chosen, the sampling rate of the "
                          "minuend should be"
                          " equal to the sampling rate of the subtrahend.")
        form.addParam('inputVolume2', PointerParam, pointerClass="Volume",
                         condition=('mapOrModel==%d ' % 0),
                         important=True, allowsNull=True,
                         label='Map to subtract (subtrahend)',
                         help="Map that has to be subtracted from the minuend 3D map.")
        form.addParam('resolution', FloatParam,
                         condition=('mapOrModel==%d' % 1),
                         label='Map resolution (A):',
                         help=" The atomic structure file wil be used to create a 3D map. "
                              " Each atom is described as a 3D Gaussian distribution of "
                              " width proportional to the resolution and amplitude proportional "
                              " to the atomic number. We recommend half the value of the"
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
                              "of the atomic structure. These removed residues "
                              "might help you to establish a control of "
                              "appropriate levels of map density.\n"
                              "In order to better visualize the area of removed "
                              "residues, 10 residues will be highligthed before "
                              "and after the first and the last residues selected,"
                              " respectively.\n")
        form.addParam('inputStructureChain', StringParam,
                         condition=(('mapOrModel==%d and '
                                    'removeResidues==True and '
                                    'selectChain==False') % 1),
                         label="Chain ", important=True,
                         help="Select a particular chain of the atomic "
                              "structure.")
        form.addParam('residuesToRemove', StringParam,
                         condition=(('mapOrModel==%d and '
                                    'removeResidues==True') % 1),
                         label="Residues to remove ", default=None,
                         important=True,
                         help="Select the first and last residues of the selected atomic "
                              "structure chain that you want to remove. \n(Use Ctrl for multiple selection)")
        form.addParam('applySymmetry', BooleanParam,
                      condition=('mapOrModel==%d' % 1),
                      label="Apply symmetry to the atomic structure:",
                      default=False,
                      help="'Symmetry = Yes' indicates that symmetry will be applied. "
                           "This option is recommended if the atomic structure "
                           "corresponds to the asymmetrical unit and you want to "
                           "regenerate the structure of the whole map.\n'Symmetry = "
                           "No' indicates that symmetry will not be applied, and then  "
                           "the map derived from the atomic structure involves only the "
                           "atomic structure provided as input.\n")
        form.addParam('symmetryGroup', EnumParam,
                      condition=(('mapOrModel==%d and applySymmetry==True') % 1),
                      choices=CHIMERA_LIST,
                      default=CHIMERA_I222r,
                      important=True,
                      label="Symmetry",
                      help="https://scipion-em.github.io/docs/docs/developer/symmetries\n"
                           "Symmetry for a description of the symmetry groups "
                           "format in CHIMERA.\n"
                           "If no symmetry is present, use _c1_."
                           'More information: \n'
                           'https://www.cgl.ucsf.edu/chimerax/docs/user/commands/sym.html'
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
                      help="This value allows to generate copies with centers "
                           "within a certain range "
                           "of distance of the center of the original molecule"
                           " model. A models's center "
                           "is defined as the center of its bounding box.")
        form.addParam('inputPdbFiles', MultiPointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      label='Other atomic structures',
                      help="In case you need to load more PDBx/mmCIF files, "
                           "you can load them here. These structures will NOT "
                           "be used "
                           "to create the subtrahend.")
        form.addParam('selectAreaMap', BooleanParam,
                      condition=('mapOrModel==%d' % 1),
                      label="Map fraction around the atomic structure?",
                      default=False,
                      help="Select 'Yes' if you want to limit the map"
                           " to a certain radius around the atomic structure.\n")
        form.addParam('radius', IntParam, default=15,
                      condition=(('mapOrModel==%d' % 1) and
                                 'selectAreaMap==True'),
                      label="Atom radius (A)",
                      help="Set the radius (Angstroms) to select values "
                           "of grid points "
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
                      condition='filterToApplyToDiffMap==%d' % 0,
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
        form.addLine(''' 
                    # vol subtract #%d #%d modelId #%d minRms true onGrid #%d
                    # (If you want to use another level the above
                    #  command recalculates the difference)
                    scipionwrite model #n [prefix stringAddedToFilename]
                    scipionss
                    scipionrs
                    scipioncombine #n1,n2,... [modelid n]
                    Type 'help command' in chimera command line for details 
                    (command is the command name)''')


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
        f.write("from chimerax.core.commands import run\n")

        # building coordinate axes
        dim = self.vol.getDim()[0]
        sampling = self.vol.getSamplingRate()
        bildFileName = os.path.abspath(self._getTmpPath("axis_input.bild"))
        Chimera.createCoordinateAxisFile(dim, bildFileName=bildFileName,
                                         sampling=sampling)
        # origin coordinates
        modelId = 1 # axis
        f.write("run(session, 'open %s')\n" % (bildFileName))
        f.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates
        # input volume
        modelMapM = modelId + 1 # 2, Minuend, result = minuend − subtrahend
        f.write("run(session,'open %s')\n" % self.fnVolName)
        # step = 1 -> no  binning
        f.write("run(session,'volume #%d style surface voxelSize %f')\n"
                % (modelMapM, sampling))
        x, y, z = self.vol.getShiftsFromOrigin()
        f.write("run(session,'volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                % (modelMapM, x, y, z))
        modelMapS =7
        modelAtomStruct = 3
        modelAtomStructChain = 4
        modelAtomStructChainSym = 5
        modelIdZone = 6
        modelMapDiff = 8
        modelMapDiffFil = 9

        if self.mapOrModel == 0:  # subtrahend is a 3D Map
            # input map
            # with its origin coordinates
            # modelMapS = modelMapM + 1  # 3 Subtrahend
            f.write("run(session,'open %s')\n" %
                    (self.subVolName))
            f.write("run(session, 'rename #3 id #%d')\n" % modelMapS)
            f.write("run(session,'volume #%d style surface voxelSize %f step 1')\n"
                    % (modelMapS, sampling))
            x, y, z = self.subVol.getShiftsFromOrigin()
            f.write("run(session,'volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                        % (modelMapS, x, y, z))
            if self.subtractOrMask == 1 and self.level.get() is not None:
                f.write("run(session,'volume #%d level %f')\n" %
                        (modelMapS, self.level))
        else:  # subtrahend is an atomic structure
            f.write("run(session,'open %s')\n" % self.atomStructName)
            # input atomic structure
            if self.selectChain == True:
                # model and chain selected
                if self.selectStructureChain.get() is not None:
                    chain = self.selectStructureChain.get()
                    self.selectedModel = chain.split(',')[0].split(':')[1].strip()
                    #TODO: Study problems with multimodels
                    if int(self.selectedModel) != 0:
                        modelId = int(modelAtomStruct +
                                                     int(self.selectedModel))

                        f.write("run(session, 'rename #%d id #%d')\n" %
                                (modelId, modelAtomStruct))
                    self.selectedChain = \
                        chain.split(',')[1].split(':')[1].strip().split('"')[1]
                    print("Selected chain: %s from model: %s from structure: %s" \
                        % (self.selectedChain, self.selectedModel,
                            os.path.basename(self.atomStructName)))
                    f.write("run(session,'sel #%d/%s')\n"
                            % (modelAtomStruct, self.selectedChain))
                    tmpPath = os.path.abspath(self._getTmpPath('chain.cif'))

                    f.write("run(session,"
                            "'save %s format mmcif models #%d relModel #%d selectedOnly true')\n"
                            % (tmpPath, modelAtomStruct, modelId))
                    f.write("run(session,'open %s')\n" % tmpPath)
                    f.write("run(session,'scipionwrite #%d prefix chain_%s_ ')\n"
                            % (modelAtomStructChain, self.selectedChain))
                    if self.selectAreaMap == True:  # mask the minuend using the atomic structure
                        if self.applySymmetry == True and self.symmetryGroup.get() is not None:
                            sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                            modelId = modelAtomStructChain #4
                            self.symMethod(f, modelId, sym, self.symmetryOrder)

                            f.write("run(session,'volume zone #%d nearAtoms #%d "
                                    "range %d newMap true modelId #%d')\n"
                                    % (modelMapM, modelAtomStructChainSym,
                                       self.radius, modelIdZone))
                            if not self.removeResidues:
                                f.write("run(session,'scipionwrite #%d prefix sym_  ')\n"
                                        % modelAtomStructChainSym)

                            f.write("run(session,'close #%d')\n" % (modelAtomStructChainSym))
                        else:
                            f.write("run(session,'volume zone #%d nearAtoms #%d "
                                    "range %d newMap true modelId #%d')\n"
                                    % (modelMapM, modelAtomStructChain,
                                    self.radius, modelIdZone))

                        f.write("run(session,'scipionwrite #%d prefix zone_  ')\n" % modelIdZone)

                    if self.removeResidues == True:
                        idxRemove = self.getIdxRemoveResidues()
                        if idxRemove is not None:
                            self.firstResidue = idxRemove[0]
                            self.lastResidue = idxRemove[1]
                            f.write("run(session,'sel #%d/%s:%d-%d')\n"
                                    % (modelAtomStructChain, self.selectedChain,
                                       int(self.firstResidue), int(self.lastResidue)))
                            f.write("run(session,'del sel')\n")
                            f.write("run(session,'sel #%d/%s:%d-%d')\n" %
                                    (modelAtomStructChain, self.selectedChain,
                                     int(self.firstResidue) - 10,
                                     int(self.lastResidue) + 10))
                    if self.applySymmetry == True:
                        if self.symmetryGroup.get() is not None:
                            sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                            modelId = modelAtomStructChain
                            self.symMethod(f, modelId, sym, self.symmetryOrder, self.rangeDist)
                            f.write("run(session,'scipionwrite #%d prefix sym_  ')\n"
                                    % modelAtomStructChainSym)
                            f.write("v=run(session,'molmap #%d %0.3f gridSpacing %f')\n"
                                    % (modelAtomStructChainSym, self.resolution, sampling))
                            f.write("run(session,'rename #%d id #7' % v.id[0])\n") ## #7 is modelMapS id
                            if self.subtractOrMask == 1 and self.level.get() is not None:
                                f.write("run(session,'volume #%d level %f')\n" %
                                        (modelMapS, self.level))
                            if self.removeResidues == True:
                                idxRemove = self.getIdxRemoveResidues()
                                if idxRemove is not None:
                                    f.write("run(session,'sel #%d:%d-%d')\n" %
                                            (modelAtomStructChainSym,
                                             int(self.firstResidue) - 10,
                                             int(self.lastResidue) + 10))

                    else:
                        f.write("v=run(session,'molmap #%d %0.3f gridSpacing %f')\n"
                                % (modelAtomStructChain, self.resolution, sampling))
                        f.write("run(session,'rename #%d id #7' % v.id[0])\n") ## #7 is modelMapS id
                        if self.subtractOrMask == 1 and self.level.get() is not None:
                            f.write("run(session,'volume #%d level %f')\n" %
                                    (modelMapS, self.level))

                    f.write("run(session,'scipionwrite #%d prefix molmap_chain%s_')\n"
                            % (modelMapS, self.selectedChain))

            else:  # use whole atomic model
                f.write("run(session,'scipionwrite #%d')\n" % modelAtomStruct)
                if self.selectAreaMap == True:
                    if self.applySymmetry == True and self.symmetryGroup.get() is not None:
                        sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                        modelId = modelAtomStruct
                        self.symMethod(f, modelId, sym, self.symmetryOrder, self.rangeDist)
                        f.write("run(session, 'rename #4 id #%d')\n" % modelAtomStructChainSym)
                        f.write("run(session,'volume zone #%d nearAtoms #%d "
                                "range %d newMap true modelId #%d')\n"
                                % (modelMapM, modelAtomStructChainSym,
                                   self.radius, modelIdZone))

                        f.write("run(session,'close #%d')\n" % (modelAtomStructChainSym))

                    else:
                        f.write("run(session,'volume zone #%d nearAtoms #%d "
                                "range %d newMap true modelId #%d')\n"
                                % (modelMapM, modelAtomStruct,
                                    self.radius, modelIdZone))

                    f.write("run(session,'scipionwrite #%d prefix zone_  ')\n" % modelIdZone)

                if self.removeResidues == True:
                    idxRemove = self.getIdxRemoveResidues()
                    if (self.inputStructureChain.get() is not None and idxRemove is not None):
                        chain = self.inputStructureChain.get()
                        self.selectedModel = chain.split(',')[0].split(':')[1].strip()
                        # TODO: Study problems with multimodels
                        if int(self.selectedModel) != 0:
                            modelId = int(modelAtomStruct +
                                          int(self.selectedModel))

                            f.write("run(session, 'rename #%d id #%d')\n" %
                                    (modelId, modelAtomStruct))
                        self.selectedChain = \
                            chain.split(',')[1].split(':')[1].strip().split('"')[1]
                        print("Selected chain: %s from model: %s from structure: %s"\
                              % (self.selectedChain, self.selectedModel,
                                 os.path.basename(self.atomStructName)))
                        f.write("run(session,'sel #%d/%s')\n"
                                % (modelAtomStruct, self.selectedChain))
                        self.firstResidue = idxRemove[0]
                        self.lastResidue = idxRemove[1]
                        f.write("run(session,'sel #%d/%s:%d-%d')\n"
                                % (modelAtomStruct, self.selectedChain,
                                   int(self.firstResidue), int(self.lastResidue)))
                        f.write("run(session,'del sel')\n")
                        f.write("run(session,'scipionwrite #%d prefix mutated_ ')\n"
                                % modelAtomStruct)

                        f.write("run(session,'sel #%d/%s:%d-%d')\n" %
                                (modelAtomStruct, self.selectedChain,
                                 int(self.firstResidue) - 10,
                                 int(self.lastResidue) + 10))

                if self.applySymmetry == True:
                    if self.symmetryGroup.get() is not None:
                        sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                        modelId = modelAtomStruct
                        self.symMethod(f, modelId, sym, self.symmetryOrder, self.rangeDist)
                        f.write("run(session, 'rename #4 id #%d')\n" % modelAtomStructChainSym)
                        f.write("run(session,'scipionwrite #%d prefix sym_ ')\n"
                                % modelAtomStructChainSym)
                        idxRemove = self.getIdxRemoveResidues()
                        if (self.inputStructureChain.get() is not None and idxRemove is not None):
                            f.write("run(session, 'sel #%d/%s:%d-%d')\n" %
                                    (modelAtomStructChainSym,
                                     self.selectedChain,
                                     int(self.firstResidue) - 10,
                                     int(self.lastResidue) + 10))
                        f.write("v=run(session,'molmap #%d %0.3f gridSpacing %f')\n"
                                % (modelAtomStructChainSym, self.resolution, sampling))
                        f.write("run(session,'rename #%d id #7' % v.id[0])\n") ## #7 is modelMapS id
                        # modelMapS = modelAtomStructChainSym + 1
                        if self.subtractOrMask == 1 and self.level.get() is not None:
                            f.write("run(session,'volume #%d level %f')\n" %
                                    (modelMapS, self.level))
                else:  # no symmetry
                    f.write("v=run(session,'molmap #%d %0.3f gridSpacing %f')\n"
                            % (modelAtomStruct, self.resolution, sampling))
                    f.write("run(session,'rename #%d id #7' % v.id[0])\n") ## #7 is modelMapS id
                    if self.subtractOrMask == 1 and self.level.get() is not None:
                        f.write("run(session,'volume #%d level %f')\n" %
                                (modelMapS, self.level))
                f.write("run(session,'scipionwrite #%d prefix molmap_  ')\n" % modelMapS)

        # Generation of the differential map
        if self.selectAreaMap == True:
            modelId = modelIdZone
        else:
            modelId = modelMapM
        if self.subtractOrMask == 0:
            f.write("run(session, 'volume subtract #%d #%d modelId #%d "
                    "minRms true onGrid #%d')\n" %
                    (modelId, modelMapS, modelMapDiff, modelId))
        else:
            f.write("run(session, 'volume mask #%d surfaces #%d invertMask "
                    "true modelId #%d')\n" %
                    (modelId, modelMapS, modelMapDiff))
        f.write("run(session,'scipionwrite #%d prefix difference_')\n" % modelMapDiff)

        # Generation of the filtered map
        if self.filterToApplyToDiffMap.get() == 0:
            f.write("run(session,'volume gaussian #%d sd %0.3f modelId %#d')\n"
                    % (modelMapDiff, self.widthFilter.get(), modelMapDiffFil))
        else:
            f.write("run(session,'volume laplacian #%d')\n" % modelMapDiff)

        f.write("run(session,'scipionwrite #%d prefix filtered_')\n" % modelMapDiffFil)
        if self.inputPdbFiles is not None:  # Other atomic models different
                                            # from the subtrahend
            for atomStruct in self.inputPdbFiles:
                f.write("run(session,'open %s')\n" %
                        os.path.abspath(atomStruct.get().getFileName()))
        # Finally save session
        f.write("run(session,'scipionss')\n")

        # run the script:
        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
            args = " --nogui --script " + \
                   os.path.abspath(self._getTmpPath(chimeraScriptFileName))
        else:
            args = " --script " + \
                   os.path.abspath(self._getTmpPath(chimeraScriptFileName))
        f.close()

        self._log.info('Launching: ' + Plugin.getProgram() + ' ' + args)

        # run in the background
        cwd = os.path.abspath(self._getExtraPath())
        Plugin.runChimeraProgram(Plugin.getProgram(), args, cwd=cwd, extraEnv=getEnvDictionary(self))

    def createOutput(self):
        # Check vol and pdb files
        directory = self._getExtraPath()
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
                kwargs = {keyword: vol}
                self._defineOutputs(**kwargs)

            elif filename.endswith(".pdb") or filename.endswith(".cif"):
                path = os.path.join(directory, filename)
                pdb = AtomStruct()
                pdb.setFileName(path)
                keyword = filename.split(".pdb")[0].replace(".","_")
                kwargs = {keyword: pdb}
                self._defineOutputs(**kwargs)

    def symMethod(self, f, modelId, sym, order=None, range=None):
        if sym == "Cn" and order != 1:
            f.write("run(session,'sym #%d C%d copies t')\n"
                    % (modelId, order))
        elif sym == "Dn" and order != 1:
            f.write("run(session,'sym #%d d%d copies t')\n"
                    % (modelId, order))
        elif sym == "T222" or sym == "TZ3":
            f.write("v=run(session,'sym #%d t,%s copies t')\n"
                    % (modelId, sym[1:]))
        elif sym == "O":
            f.write("run(session,'sym #%d O copies t')\n"
                    % modelId)
        elif sym == "I222" or sym == "I222r" or sym == "In25" or \
                sym == "In25r" or sym == "I2n3" or sym == "I2n3r" or \
                sym == "I2n5" or sym == "I2n5r":
            f.write("run(session,'sym #%d i,%s copies t')\n"
                    % (modelId, sym[1:]))

        f.write("run(session,'delete #%d & #%d #>%d')\n"
                    % (int(modelId) + 1, modelId, self.rangeDist))

    def getIdxRemoveResidues(self):
        resJson = getattr(self, 'residuesToRemove').get()
        if resJson:
            idxs = json.loads(resJson)['index'].split('-')
            return idxs


    def _summary(self):
        summary = []
        if self.getOutputsSize() > 0:
            directory = self._getExtraPath()
            summary.append("Produced files:")
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".pdb"):
                    summary.append(filename)
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".mrc"):
                    summary.append(filename)
            summary.append("we have some result")
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        return summary

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

