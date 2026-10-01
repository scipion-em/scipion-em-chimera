import os
from pyworkflow.viewer import Viewer
from pwem.viewers.viewer_chimera import Chimera
from ..protocols.protocol_discrepancies import ChimeraProtDiscrepancies


class ChimeraProtDiscrepanciesViewer(Viewer):
    """ Viewer for ChimeraProtDiscrepancies protocol output. """
    _label = 'viewer discrepancies'
    _targets = [ChimeraProtDiscrepancies]

    def visualize(self, obj, **args):
        # Create Chimera command file
        fnCmd = self.protocol._getExtraPath("discrepancies_viewer_with_files.cxc")
        with open(fnCmd, 'w') as f:
            # Process protocol outputs
            for output in self.protocol._outputs:
                # If the file is an atomic structure (.cif or .pdb), open it in Chimera
                fileName = os.path.abspath(eval(f'self.protocol.{output}.getFileName()'))
                if fileName.endswith(".cif") or fileName.endswith(".pdb"):
                    f.write(f"open {fileName}\n")
                    # Apply occupancy color palette to the structure
                    f.write("color byattribute occupancy palette paegreen\n")
                    f.write("key darkgreen:low green: lightgreen: white:high\n") # This palette means RMSD, high is worse

        # Run Chimera with the generated command file
        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")