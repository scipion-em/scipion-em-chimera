# **************************************************************************
# *
# * Authors:  Roberto Marabini (roberto@cnb.csic.es), May 2013
# *           Marta Martinez (mmmtnez@cnb.csic.es)
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

from pwem.convert import Ccp4Header
from pwem.emlib.image import ImageHandler
from pwem.objects import Volume
from pwem.objects import Transform
import pyworkflow.viewer as pwviewer

from ..protocols import ChimeraSubtractionMaps
from ..protocols.protocol_fit import ChimeraProtRigidFit
from ..protocols.protocol_operate import ChimeraProtOperate
from ..protocols.protocol_restore import ChimeraProtRestore
from ..protocols.protocol_modeller_search import ChimeraModelFromTemplate
from ..protocols.protocol_alphafold import ChimeraImportAtomStructAlphafold

from pwem.viewers.viewer_chimera import (Chimera,
                                         sessionFile)
from pyworkflow.viewer import DESKTOP_TKINTER, Viewer
from chimera.objects import PAE

class ChimeraViewerBase(Viewer):
    """ Visualize the output of protocols protocol_fit and protocol_operate """
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **args):
        # The input map or pdb may be a parameter from the protocol
        # or from the parent protocol.
        dim = 150.
        sampling = 1.
        _inputVol = None
        directory = self.protocol._getExtraPath()
        try:
            try:
                if self.protocol.inputVolume.get() is not None:
                    _inputVol = self.protocol.inputVolume.get()
                elif self.protocol.pdbFileToBeRefined.get().getVolume() is not None:
                    _inputVol = self.protocol.pdbFileToBeRefined.get().getVolume()
                elif self.protocol.inputVolumes[0] is not None:
                    _inputVol = self.protocol.inputVolumes[0].get()
            except:
                output3DMapList = []
                for filename in sorted(os.listdir(directory)):
                    if filename.endswith(".mrc"):
                        output3DMapList.append(filename.split('.')[0])
                        output3DMap = str(output3DMapList[0])
                        if len(output3DMap) > 0:
                            _inputVol = self.protocol.output3DMap
        except:
            # TODO: I do not know if we still need this part
            # Remark that inputProtocol does not longer exist, it has been replaced by inputProtocolDict
            # Compare with the previous code, specially the alternative directory
            for item in list(self.protocol.inputProtocolDict().values()):
                if item.hasAttribute('inputVolume') and item.inputVolume.get() is not None:
                    _inputVol = item.inputVolume.get()
                    break
                elif item.hasAttribute('pdbFileToBeRefined') and \
                    item.pdbFileToBeRefined.get().getVolume() is not None:
                    _inputVol = item.pdbFileToBeRefined.get().getVolume()
                    break
                # directory = item._getExtraPath()

        if _inputVol is not None:
            dim = _inputVol.getDim()[0]
            sampling = _inputVol.getSamplingRate()

        bildFileName = self.protocol._getExtraPath("axis_output.bild")
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=bildFileName,
                                         sampling=sampling)
        fnCmd = self.protocol._getExtraPath("chimera_output.cxc")
        f = open(fnCmd, 'w')
        # change to workingDir
        # If we do not use cd and the project name has an space
        # the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())
        f.write("open %s\n" % bildFileName)
        f.write("cofr 0,0,0\n")  # set center of coordinates
        inputVolFileName = ''
        counter = 2
        if _inputVol is not None:
            # In case we have PDBs only, _inputVol is None:
            self.visInputVolume(f, _inputVol, counter)
        else:
            counter = 1

        if (self.protocol.hasAttribute("inputVolume2") and\
                self.protocol.inputVolume2.get() is not None):
            counter += 1
            _inputVol2 = self.protocol.inputVolume2.get()
            self.visInputVolume(f, _inputVol2, counter)

        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".mrc") and filename != inputVolFileName:
                counter += 1
                volFileName = os.path.join(directory, filename)
                vol = Volume()
                vol.setFileName(volFileName)

                # fix mrc header
                ccp4header = Ccp4Header(volFileName, readHeader=True)
                sampling = ccp4header.computeSampling()
                origin = Transform()
                shifts = ccp4header.getOrigin()
                origin.setShiftsTuple(shifts)
                f.write("open %s\n" % volFileName)
                f.write("volume #%d style surface voxelSize %f\n"
                        "volume #%d origin %0.2f,%0.2f,%0.2f\n"
                        % (counter, sampling, counter, shifts[0], shifts[1], shifts[2]))
                f.write("volume #%d level %0.3f\n"
                        % (counter, 0.001))

        for filename in os.listdir(directory):
            if filename.endswith(".pdb") or filename.endswith(".cif"):
                if not (filename.startswith("Atom_struct_out_") or
                        filename.startswith("tmp_")):
                    path = os.path.join(directory, filename)
                    f.write("open %s\n" % path)
                    
        f.write("view\n")
        f.close()

        # run in the background
        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")
        return []

    def visInputVolume(self, f, vol, counter):
        inputVolFileName = ImageHandler.removeFileType(vol.getFileName())
        f.write("open %s\n" % inputVolFileName)
        if vol.hasOrigin():
            x, y, z = vol.getOrigin().getShifts()
        else:
            x, y, z = vol.getOrigin(force=True).getShifts()
        f.write("volume #%d style surface voxelSize %f\n"
                "volume #%d origin %0.2f,%0.2f,%0.2f\n"
                % (counter, vol.getSamplingRate(), counter, x, y, z))

class ChimeraRestoreViewer(Viewer):
    """ Visualize the output of protocols protocol_fit and protocol_operate """
    _label = 'viewer restore'
    _targets = [ChimeraProtRestore]

    def _visualize(self, obj, **args):
        fnCmd = self.protocol._getExtraPath("chimera_restore_session.cxc")
        f = open(fnCmd, 'w')
        # change to workingDir
        # If we do not use cd and the project name has an space
        # the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())
        path1 = os.path.join(self.protocol._getExtraPath(), sessionFile)
        if os.path.exists(path1):
            # restored SESSION
            path = path1
        else:
            # SESSION from inputProtocol
            path2 = os.path.join(
                self.protocol.inputProtocol.get()._getExtraPath(), sessionFile)
            path = path2
        f.write("open %s\n" % path)

        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")
        return []


class ChimeraProtRigidFitViewer(ChimeraViewerBase):
    _label = 'viewer fit'
    _targets = [ChimeraProtRigidFit]


class ChimeraProtOperateViewer(ChimeraViewerBase):
    _label = 'viewer operate'
    _targets = [ChimeraProtOperate]


class ChimeraModelFromTemplateViewer(ChimeraViewerBase):
    _label = 'viewer model from template'
    _targets = [ChimeraModelFromTemplate]

class ChimeraSubtractionMapsViewer(ChimeraViewerBase):
    _label = 'viewer subtract maps'
    _targets = [ChimeraSubtractionMaps]

class PAEViewer(Viewer):
    """Visulize Alphafold PAE objects
    object defined at from chimera.objects import PAE
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [PAE]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)

    def visualize(self, obj, **kwargs):
        paeObject = obj
        paeObject.read()
        paeMatrix = paeObject.getMatrix()
        import matplotlib.pyplot as plt
        if paeMatrix is None:
            raise("matrix is None. I cannot plot an empty matrix")
        else:
            fig2 = plt.figure(2)
            plt.title('Predicted aligned error\n Expected position error (Ångströms)')
            plt.xlabel('Scored residue')
            plt.ylabel('Aligned residue')
            plt.imshow(paeMatrix, interpolation='none', cmap='Greens_r')
            plt.colorbar()
        plt.show()

class ChimeraAlphafoldViewer(Viewer):
    _label = 'viewer alphafold'
    _targets = [ChimeraImportAtomStructAlphafold]
    def _visualize(self, obj, **args):
        # create axis file
        models = 1
        dim = 150
        sampling = 1.
        extraFileName = os.path.abspath(self.protocol._getExtraPath("axis_input.bild"))
        Chimera.createCoordinateAxisFile(dim,
                                         bildFileName=extraFileName,
                                         sampling=sampling)

        fnCmd = self.protocol._getExtraPath("chimera_alphafold.cxc")
        f = open(fnCmd, 'w')
        f.write("open %s\n" % extraFileName)
        models +=1
        f.write("cofr 0,0,0\n")  # set center of coordinates
        # change to workingDir
        # If we do not use cd and the project name has an space
        # the protocol fails even if we pass absolute paths
        f.write('cd %s\n' % os.getcwd())

        # get path to atomstructs
        for output in self.protocol._outputs:
            # if the file is an atomic struct show it in chimera
            fileName = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
            if fileName.endswith(".cif") or fileName.endswith(".pdb"):
                f.write("open %s\n" % fileName)
                models +=1
        # if exists upload other results files 
        # model_?_unrelaxed.pdb
        #pattern = self.protocol._getExtraPath("results/model_?_unrelaxed.pdb")
        #from glob import glob
        #for model in glob(pattern):
        #    f.write("open %s\n" % model)
        #    f.write(f"hide #{models} models\n")
        #    models +=1
        # set alphafold colormap
        f.write("color bfactor palette alphafold\n")
        f.write("key red:low orange: yellow: cornflowerblue: blue:high\n")
        f.close()
        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")
        # plot coverage
        if self.protocol.source == self.protocol.IMPORT_REMOTE_ALPHAFOLD and \
           (self.protocol.colabID.get() == self.protocol.CHIMERA21 or
            self.protocol.colabID.get() == self.protocol.TEST):
            database_names=[]
            if self.protocol.colabID.get() == self.protocol.CHIMERA21:
                database_names=["sequence_1_mgnify", 
                                "sequence_1_smallbfd", 
                                "sequence_1_uniref90"]
                paeFile = "best_model_pae.json"
            elif self.protocol.colabID.get() == self.protocol.TEST:
                database_names=["sequence_1_mgnify", 
                                "sequence_1_smallbfd", 
                                "sequence_1_uniref90"]
                paeFile = "best_model_pae.json"
            self.plot_alignment_coverage(database_names=database_names,
                                         paeFile = paeFile)
        elif self.protocol.colabID.get() == self.protocol.PHENIX:
            for output in self.protocol._outputs:
                # if the file is an atomic struct show it in chimera
                fileName = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
                if fileName.endswith(".jsn"):
                    paeFile = os.path.basename(fileName)
            self.plot_alignment_coverage(database_names=None,
                                         paeFile = paeFile)
        return []

    def plot_alignment_coverage(self, database_names=["mgnify", "smallbfd", "uniref90"], paeFile = None):
        def read_alignments(database_names, output_dir):
            alignments, deletions = [], []
            from os import path
            for name in database_names:
                apath = path.join(output_dir, name + '_alignment')
                dpath = path.join(output_dir, name + '_deletions')
                if not path.exists(apath) or not path.exists(dpath):
                    return [],[]
                with open(apath, 'r') as f:
                    seqs = [line.rstrip() for line in f.readlines()]
                    alignments.append(seqs)
                with open(dpath, 'r') as f:
                    dcounts = [[int(value) for value in line.split(',')] for line in f.readlines()]
                    deletions.append(dcounts)
            return alignments, deletions

        def _plot_alignment_coverage(alignments, paeMatrix):
            import matplotlib.pyplot as plt
            if alignments is not None:
                counts = alignment_coverage(alignments)
                if counts is None:
                    return
                fig1 = plt.figure(1, figsize=(12, 3))
                plt.title('Number of Aligned Sequences with no Gap for each Residue Position')
                x = range(1, len(counts)+1) # Start residue numbers at 1, not 0.
                plt.plot(x, counts, color='black')
                plt.xlabel('Residue number')
                plt.ylabel('Coverage')

    
            if paeMatrix is None:
                raise("matrix is None. I cannot plot an empty matrix")
            else:
                fig2 = plt.figure(2)
                plt.title('Predicted aligned error\n Expected position error (Ångströms)')
                plt.xlabel('Scored residue')
                plt.ylabel('Aligned residue')
                plt.imshow(paeMatrix, interpolation='none', cmap='Greens_r')
                plt.colorbar()
            plt.show()


        def alignment_coverage(alignments):
            from numpy import zeros, int32
            counts = None
            for alignment in alignments:
                for line in alignment:
                    if counts is None:
                        counts = zeros((len(line),), int32)
                    for i,c in enumerate(line):
                        if c != '-':
                            counts[i] += 1
            return counts

        output_dir = self.protocol._getExtraPath('results')
        if database_names is not None:
            alignments, deletions = read_alignments(database_names, output_dir)
        else:
            alignments = None
        paeFile = os.path.abspath(self.protocol._getExtraPath(os.path.join('results', paeFile)))
        paeObject = PAE(filename=paeFile)
        paeObject.read()
        paeMatrix = paeObject.getMatrix()
        _plot_alignment_coverage(alignments, paeMatrix)
