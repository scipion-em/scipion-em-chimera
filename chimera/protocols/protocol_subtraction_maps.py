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
from .protocol_base import createScriptFile, ChimeraProtBase
from pwem.viewers.viewer_chimera import (Chimera,
                                         sessionFile,
                                         chimeraMapTemplateFileName,
                                         chimeraScriptFileName,
                                         chimeraPdbTemplateFileName)
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from .. import Plugin
from pyworkflow.utils.properties import Message

class ChimeraSubtractionMaps(EMProtocol):
    """Protocol to subtract two volumes.
        One of these volumes can be derived from an atomic structure.
        Execute command *scipionwrite model #n [refmodel #p] [prefix stringAddedToFilename]*
        from command line in order to transfer the generated maps and models to scipion.
        In addition to maps and models that the protocol saves by default,
        the user can generate and save some others"""
    _label = 'map subtraction'
    _program = ""
    _program = ""
    _version = VERSION_3_0

    PROTOCOL_OPTIONS = ['Subtraction', 'Mask']
    MAP_OPTIONS = ['3D map', 'atomic structure']
    CHIMERA_FILTERS = ['Gaussian', 'Fourier Transform']

    subtractionString = """
from VolumeStatistics import mean_sd_rms
from chimera import openModels
from chimera import replyobj

from VolumeData import Array_Grid_Data
from VolumeViewer.volume import volume_from_grid_data
import numpy
from numpy import greater_equal, multiply, dot as inner_product
from numpy import array, ravel
# get volume from model id
def subtraction(minuendId, subtrahendId, outModelId=-1, subtractOrMask=0):
    ''' subtract or mask two volumes after adjust their respective ranges.
    
        The mask part is close to the chimera command
        vop zone invert. 
        If subtractOrMask=1 the program function uses vop subtract (chimera) 
        Else:
        1) A mask is computed using the volume with modelid=subtrahendId 
           and the countour level value.
        2) All voxels in minuendId are set to 0 for those
           voxels inside the mask.
        3) Result is shown in chimera
        4) Volumes are assume to have the same sampling rate and dimensions
        Usage Example:
        from chimera import runCommand
        runCommand('open /home/roberto/Downloads/Vols/emd_21375_crop_ref.mrc')
        runCommand('open /home/roberto/Downloads/Vols/i2pc_Level0_226_crop_ref.mrc')
        subtraction(0, 1, outModelId=6, subtractOrMask=0)
        
        Note this function is milar to chimera's vop zone invert
        but can be used not only with PDBs but with 3D maps as subtrahend
    '''
    if subtractOrMask==0:
        command = "vop subtract #%d #%d modelId #%d minRMS true onGrid #%d" % (minuendId, subtrahendId, outModelId, minuendId)
        runCommand(command)
        return
    
    # get models from Ids
    minuendModel = openModels.list(id=minuendId)[0]  #submodel if needed
    subtrahendModel = openModels.list(id=subtrahendId)[0]  #submodel if needed
    
    # Get contour level from model
    contourLevel = subtrahendModel.surface_levels[0]
        
    # get  matrix with voxel values
    minuendMatrix = minuendModel.full_matrix()
    subtrahendMatrix = subtrahendModel.full_matrix()
    
    # check if sampling and size is the same
    compatible1 = abs (minuendModel.data.step[0] - minuendModel.data.step[0]) < 0.01
    s1 = minuendMatrix.shape
    s2 = subtrahendMatrix.shape
    compatible2 = (s1[0]==s2[0]) & (s1[1]==s2[1]) & (s1[1]==s2[1])

    if not (compatible1 & compatible2):
        replyobj.status("Both volumes have incompatible size or sampling, using vop subtract")
        command = "vop subtract #%d #%d modelId #%d minRMS true onGrid #%d" % (minuendId, subtrahendId, outModelId, minuendId)
        runCommand(command)
        return

    # test data, comment next two lines to operate
    # with chimera volumes
    # minuendMatrix = array([[1., 2., 3.], [4., 5., 6.]])
    # subtrahendMatrix = minuendMatrix * 2 +1
    # contourLevel=6
    
    
    # create mask
    # keep only values above threshold contourLevel
    mask = subtrahendMatrix > contourLevel  # matrix with true and false
                                            # values less than counter
                                            # are set to true
    # set values to half the contourLevel                                           
    minuendMatrix[mask] = 0.
    
    #innerProduct = inner_product(shifted_matrix_subtrahend[mask],    shifted_matrix_minuend[mask])
    #normalization =inner_product(shifted_matrix_subtrahend[mask], shifted_matrix_subtrahend[mask])
    #if normalization == 0:
    #   f = 1
    #else:
    #   f = innerProduct/normalization
    #new_matrix_subtrahend =  minuendMatrix - (shifted_matrix_subtrahend * f + minuendMean )
    # attach matrix to grid
    
    g0 = Array_Grid_Data(minuendMatrix, minuendModel.data.origin, 
                         minuendModel.data.step, minuendModel.data.cell_angles)
    # create new volume
    if outModelId == -1:
        differenceVolume = volume_from_grid_data(g0)
    else:
        differenceVolume = volume_from_grid_data(g0, model_id=outModelId)
"""

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form, doHelp=False):

        form.addSection(label='Input')
        form.addParam('inputVolume', PointerParam, pointerClass="Volume",
                      label='Input 3D Map',
                      important=True,
                      help="difference 3D map = minuend − subtrahend"
                           "input here the minuend")
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
                      help='result = minuend − subtrahend. Calculation are made for those\n'
                           'voxels inside a region created by this contour level\n'
                           'empty -> chimera computes the level')
        form.addParam('mapOrModel', EnumParam,
                     choices=self.MAP_OPTIONS,
                     display=EnumParam.DISPLAY_HLIST,
                     default=0, label='Subtraction/Mask of',
                     help="difference 3D Map = minuend − subtrahend. "
                          "Subtrahend 3D map may be provided by the user (choose '3D map') "
                          "or created from an atomic coordinates file (choose 'atomic structure'"
                          " If 3D Map is chosen, the sampling rate of the minuend should be"
                          " equal to the sampling rate of the subtrahend.")
        form.addParam('inputVolume2', PointerParam, pointerClass="Volume",
                         condition=('mapOrModel==%d ' % 0),
                         important=True, allowsNull=True,
                         label='Map to subtract (subtrahend)',
                         help="Map that has to be subtracted from the minuend 3D map")
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
                              "In order to better visualize the area of removed"
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
                      default=CHIMERA_I222r,
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
                    vop subtract #%d #%d modelId #%d minRMS true onGrid #%d
                    (If you want to use another level the above
                     command recalculates the difference)
                    scipionwrite model #n [refmodel #p] [prefix stringAddedToFilename]
                    scipionss
                    scipionrs
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
        # origin coordinates
        modelId = 0 # axis
        f.write("runCommand('open %s')\n" % (bildFileName))
        f.write("runCommand('cofr 0,0,0')\n")  # set center of coordinates
        # input volume
        modelMapM = modelId + 1 # 1, Minuend, result = minuend − subtrahend
        f.write("runCommand('open %s')\n" % self.fnVolName)
        # step = 1 -> no  binning
        f.write("runCommand('volume #%d style surface voxelSize %f')\n"
                % (modelMapM, sampling))
        x, y, z = self.vol.getShiftsFromOrigin()
        f.write("runCommand('volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                % (modelMapM, x, y, z))

        if self.mapOrModel == 0:  # subtrahend is a 3D Map
            # input map
            # with its origin coordinates
            modelMapS = modelMapM + 1  # 2 Subtrahend
            f.write("runCommand('open %s')\n" %
                    (self.subVolName))
            f.write("runCommand('volume #%d style surface voxelSize %f step 1')\n"
                    % (modelMapS, sampling))
            x, y, z = self.subVol.getShiftsFromOrigin()
            f.write("runCommand('volume #%d origin %0.2f,%0.2f,%0.2f')\n"
                        % (modelMapS, x, y, z))
            if self.subtractOrMask == 1 and self.level.get() is not None:
                f.write("runCommand('volume #%d level %f')\n" %
                        (modelMapS, self.level))
        else:  # subtrahend is an atomic structure
            f.write("runCommand('open %s')\n" % self.atomStructName)
            # input atomic structure
            modelAtomStruct = modelMapM + 1
            if self.selectChain == True:
                # model and chain selected
                if self.selectStructureChain.get() is not None:
                    chain = self.selectStructureChain.get()
                    self.selectedModel = chain.split(',')[0].split(':')[1].strip()
                    modelAtomStruct = int(modelAtomStruct +
                                                 int(self.selectedModel))
                    self.selectedChain = \
                        chain.split(',')[1].split(':')[1].strip().split('"')[1]
                    print("Selected chain: %s from model: %s from structure: %s" \
                        % (self.selectedChain, self.selectedModel,
                            os.path.basename(self.atomStructName)))
                    f.write("runCommand('sel #%d:.%s')\n"
                            % (modelAtomStruct, self.selectedChain))
                    tmpPath = self._getTmpPath('chain.pdb')
                    f.write("runCommand('write format pdb selected relative %d #%d %s')\n"
                            % (modelId, modelAtomStruct, tmpPath))
                    f.write("runCommand('open %s')\n" % tmpPath)
                    modelAtomStructChain = modelAtomStruct + 1
                    f.write("runCommand('scipionwrite model #%d refmodel #%d "
                            "prefix chain_%s_  savesession 0')\n"
                            % (modelAtomStructChain, modelMapM,
                               self.selectedChain))
                    if self.selectAreaMap == True:  # mask the minuend using the atomic structure
                        if self.applySymmetry == True and self.symmetryGroup.get() is not None:
                            sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                            modelId = modelAtomStructChain
                            self.symMethod(f, modelId, sym, self.symmetryOrder, self.rangeDist)

                            modelAtomStructChainSym = modelAtomStructChain + 2
                            # create a new model with the result of the symmetrization
                            f.write("runCommand('combine #%d- modelId #%d')\n"
                                    % (modelAtomStructChain, modelAtomStructChainSym))
                            modelIdZone = modelAtomStructChainSym + 1
                            f.write("runCommand('vop zone #%d #%d %d modelId #%d')\n"
                                    % (modelMapM, modelAtomStructChainSym,
                                       self.radius, modelIdZone))

                            f.write("runCommand('close #%d')\n" % (modelAtomStructChain + 1))
                            f.write("runCommand('close #%d')\n" % (modelAtomStructChainSym))
                        else:
                            modelIdZone = modelAtomStructChain + 1
                            f.write("runCommand('vop zone #%d #%d %d modelId #%d')\n"
                                    % (modelMapM, modelAtomStructChain,
                                    self.radius, modelIdZone))

                        f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                                "prefix zone_  savesession 0')\n" % (modelIdZone, modelMapM))
                        modelMapS = modelIdZone + 1
                    else:  # do not mask the minuend using the atomic structure
                        if self.applySymmetry == True:
                            modelMapS = modelAtomStructChain + 3
                        else:
                            modelMapS = modelAtomStructChain + 1

                    if self.removeResidues == True:
                        if (self.firstResidueToRemove.get() is not None and
                                self.lastResidueToRemove.get() is not None):
                            self.firstResidue = self.firstResidueToRemove.get().\
                            split(":")[1].split(",")[0].strip()
                            self.lastResidue = self.lastResidueToRemove.get(). \
                                split(":")[1].split(",")[0].strip()
                            f.write("runCommand('select #%d:%d-%d.%s')\n"
                                    % (modelAtomStructChain,
                                       int(self.firstResidue), int(self.lastResidue),
                                       self.selectedChain))
                            f.write("runCommand('del sel')\n")
                            f.write("runCommand('select #%d:%d-%d.%s')\n" %
                                    (modelAtomStructChain, int(self.firstResidue) - 10,
                                     int(self.lastResidue) + 10, self.selectedChain))
                    if self.applySymmetry == True:
                        if self.symmetryGroup.get() is not None:
                            sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                            modelId = modelAtomStructChain
                            self.symMethod(f, modelId, sym, self.symmetryOrder, self.rangeDist)
                            modelAtomStructChainSym = modelAtomStructChain + 2
                            f.write("runCommand('combine #%d- modelId #%d')\n"
                                    % (modelAtomStructChain, modelAtomStructChainSym))
                            f.write("runCommand('scipionwrite model #%d refmodel #%d "
                                    "prefix sym_  savesession 0')\n"
                                    % (modelAtomStructChainSym, modelMapM))
                            f.write("runCommand("
                                    "'molmap #%d %0.3f gridSpacing %f modelId #%d')\n"
                                    % (modelAtomStructChainSym, self.resolution, sampling,
                                       modelMapS))
                            if self.subtractOrMask == 1 and self.level.get() is not None:
                                f.write("runCommand('volume #%d level %f')\n" %
                                        (modelMapS, self.level))
                            if self.removeResidues == True:
                                if (self.firstResidueToRemove.get() is not None and
                                        self.lastResidueToRemove.get() is not None):
                                    f.write("runCommand('select #%d:%d-%d')\n" %
                                            (modelAtomStructChainSym,
                                             int(self.firstResidue) - 10,
                                             int(self.lastResidue) + 10))

                    else:
                        f.write("runCommand("
                                "'molmap #%d %0.3f gridSpacing %f modelId #%d')\n"
                                % (modelAtomStructChain, self.resolution, sampling,
                                   modelMapS))
                        if self.subtractOrMask == 1 and self.level.get() is not None:
                            f.write("runCommand('volume #%d level %f')\n" %
                                    (modelMapS, self.level))
                    f.write("runCommand('scipionwrite model #%d refmodel #%d "
                            "prefix molmap_chain%s_  savesession 0')\n"
                            % (modelMapS, modelMapM,
                               self.selectedChain))

            else:  # use whole atomic model
                f.write("runCommand('scipionwrite model #%d refmodel #%d  savesession 0')\n" \
                        % (modelAtomStruct, modelMapM))
                if self.selectAreaMap == True:
                    if self.applySymmetry == True and self.symmetryGroup.get() is not None:
                        sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                        modelId = modelAtomStruct
                        self.symMethod(f, modelId, sym, self.symmetryOrder, self.rangeDist)
                        modelAtomStructChainSym = modelAtomStruct + 2
                        f.write("runCommand('combine #%d- modelId #%d')\n"
                                % (modelAtomStruct, modelAtomStructChainSym))

                        modelIdZone = modelAtomStructChainSym + 1
                        f.write("runCommand('vop zone #%d #%d %d modelId #%d')\n"
                                % (modelMapM, modelAtomStructChainSym,
                                   self.radius, modelIdZone))

                        f.write("runCommand('close #%d')\n" % (modelAtomStructChainSym))

                    else:
                        modelIdZone = modelAtomStruct + 1
                        f.write("runCommand('vop zone #%d #%d %d modelId #%d')\n"
                                % (modelMapM, modelAtomStruct,
                                    self.radius, modelIdZone))

                    f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                            "prefix zone_  savesession 0')\n" % (modelIdZone,
                                                  modelMapM))
                    modelMapS = modelIdZone + 1

                else:
                    if self.applySymmetry == True:
                        modelMapS = modelAtomStruct + 3
                    else:
                        modelMapS = modelAtomStruct + 1

                if self.removeResidues == True:
                    if (self.inputStructureChain.get() is not None and
                            self.firstResidueToRemove.get() is not None and
                            self.lastResidueToRemove.get() is not None):
                        chain = self.inputStructureChain.get()
                        self.selectedModel = chain.split(',')[0].split(':')[1].strip()
                        modelAtomStruct = int(modelAtomStruct +
                                                     int(self.selectedModel))
                        self.selectedChain = \
                            chain.split(',')[1].split(':')[1].strip().split('"')[1]
                        print("Selected chain: %s from model: %s from structure: %s" \
                              % (self.selectedChain, self.selectedModel,
                                 os.path.basename(self.atomStructName)))
                        f.write("runCommand('sel #%d:.%s')\n"
                                % (modelAtomStruct, self.selectedChain))
                        self.firstResidue = self.firstResidueToRemove.get(). \
                            split(":")[1].split(",")[0].strip()
                        self.lastResidue = self.lastResidueToRemove.get(). \
                            split(":")[1].split(",")[0].strip()
                        f.write("runCommand('select #%d:%d-%d.%s')\n"
                                % (modelAtomStruct,
                                   int(self.firstResidue), int(self.lastResidue),
                                   self.selectedChain))
                        f.write("runCommand('del sel')\n")
                        f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                                "prefix mutated_  savesession 0')\n"
                                % (modelAtomStruct, modelMapM))

                        f.write("runCommand('select #%d:%d-%d.%s')\n" %
                                (modelAtomStruct, int(self.firstResidue) - 10,
                                 int(self.lastResidue) + 10, self.selectedChain))

                if self.applySymmetry == True:
                    if self.symmetryGroup.get() is not None:
                        sym = CHIMERA_SYM_NAME[self.symmetryGroup.get()]
                        modelId = modelAtomStruct
                        self.symMethod(f, modelId, sym, self.symmetryOrder, self.rangeDist)
                        modelAtomStructChainSym = modelAtomStruct + 2
                        f.write("runCommand('combine #%d- modelId #%d')\n"
                                % (modelAtomStruct, modelAtomStructChainSym))
                        f.write("runCommand('scipionwrite model #%d refmodel #%d "
                                "prefix sym_  savesession 0')\n"
                                % (modelAtomStructChainSym, modelMapM))
                        if (self.inputStructureChain.get() is not None and
                                self.firstResidueToRemove.get() is not None and
                                self.lastResidueToRemove.get() is not None):
                            f.write("runCommand('select #%d:%d-%d.%s')\n" %
                                    (modelAtomStructChainSym,
                                     int(self.firstResidue) - 10,
                                     int(self.lastResidue) + 10,
                                     self.selectedChain))
                        f.write("runCommand("
                                "'molmap #%d %0.3f gridSpacing %f modelId #%d')\n"
                                % (modelAtomStructChainSym, self.resolution, sampling,
                                   modelMapS))
                        if self.subtractOrMask == 1 and self.level.get() is not None:
                            f.write("runCommand('volume #%d level %f')\n" %
                                    (modelMapS, self.level))
                else:  # no symmetry
                    f.write("runCommand("
                            "'molmap #%d %0.3f gridSpacing %f modelId #%d')\n"
                            % (modelAtomStruct, self.resolution, sampling,
                               modelMapS))
                    if self.subtractOrMask == 1 and self.level.get() is not None:
                        f.write("runCommand('volume #%d level %f')\n" %
                                (modelMapS, self.level))
                f.write("runCommand('scipionwrite model #%d refmodel #%d "
                        "prefix molmap_  savesession 0')\n"
                        % (modelMapS, modelMapM))

        # Generation of the differential map
        f.write(self.subtractionString)

        modelMapDiff = modelMapS + 1
        if self.selectAreaMap == True:
            f.write("subtraction(%d, %d, outModelId=%d, subtractOrMask=%d)\n" %
                    (modelIdZone, modelMapS, modelMapDiff, self.subtractOrMask.get())
                    )
            #f.write("runCommand('vop subtract #%d #%d modelId #%d "
            #        "minRMS true onGrid #%d')\n"
            #        % (modelIdZone, modelMapS,
            #           modelMapDiff, modelMapM))
        else:
            f.write("subtraction(%d, %d, outModelId=%d, subtractOrMask=%r)\n" %
                    (modelMapM, modelMapS, modelMapDiff, self.subtractOrMask.get())
                    )
            #f.write("runCommand('vop subtract #%d #%d modelId #%d "
            #        "minRMS true onGrid #%d')\n"
            #        % (modelMapM, modelMapS,
            #           modelMapDiff, modelMapM))


        f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                "prefix difference_  savesession 0')\n"
                % (modelMapDiff, modelMapM))

        # Generation of the filtered map
        modelMapDiffFil = modelMapDiff + 1
        if self.filterToApplyToDiffMap.get() == 0:
            f.write("runCommand('vop gaussian #%d sd %0.3f')\n"
                    % (modelMapDiff, self.widthFilter.get()))
        else:
            f.write("runCommand('vop laplacian #%d')\n"
                    % (modelMapDiff))
        f.write("runCommand('scipionwrite model #%d refmodel #%d " \
                "prefix filtered_  savesession 0')\n"
                % (modelMapDiffFil, modelMapM))
        if self.inputPdbFiles is not None:  # Other atomic models different
                                            # from the subtrahend
            for atomStruct in self.inputPdbFiles:
                f.write("runCommand('open %s')\n" %
                        os.path.abspath(atomStruct.get().getFileName()))
        # Finally save session
        f.write("runCommand('scipionss')\n")

        # run the script:
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
        # cHB = ChimeraProtBase()
        # directory = self._getExtraPath()
        # cHB.createOutput(directory)
        """ Register outputs.
        """
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
            f.write("runCommand('sym #%d group C%d range %d')\n"
                    % (modelId, order, range))
        elif sym == "Dn" and order != 1:
            f.write("runCommand('sym #%d group d%d range %d')\n"
                    % (modelId, order, range))
        elif sym == "T222" or sym == "TZ3":
            f.write("runCommand('sym #%d group t,%s range %d')\n"
                    % (modelId, sym[1:], range))
        elif sym == "O":
            f.write("runCommand('sym #%d group O range %d')\n"
                    % modelId, range)
        elif sym == "I222" or sym == "I222r" or sym == "In25" or \
                sym == "In25r" or sym == "I2n3" or sym == "I2n3r" or \
                sym == "I2n5" or sym == "I2n5r":
            f.write("runCommand('sym #%d group i,%s range %d')\n"
                    % (modelId, sym[1:], range))

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

