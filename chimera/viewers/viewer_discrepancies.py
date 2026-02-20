import os
from pyworkflow.viewer import Viewer
from pwem.viewers.viewer_chimera import Chimera
from ..protocols.protocol_discrepancies import ChimeraProtDiscrepancies


class ChimeraProtDiscrepanciesViewer(Viewer):
    """ Viewer for ChimeraProtDiscrepancies protocol output. """
    _label = 'viewer discrepancies'
    _targets = [ChimeraProtDiscrepancies]

    def visualize(self, obj, **args):
        fnCmd = self.protocol._getExtraPath("discrepancies_viewer_with_files.cxc")
        with open(fnCmd, 'w') as f:
            outputs = self.protocol._outputs
            if not outputs:
                return

            # Reference model = first with 'ref_' or first file
            ref_file = None
            for output in outputs:
                file_path = os.path.abspath(getattr(self.protocol, output).getFileName())
                if os.path.basename(file_path).startswith("ref_"):
                    ref_file = file_path
                    break
            if not ref_file:
                ref_file = os.path.abspath(getattr(self.protocol, outputs[0]).getFileName())

            # Abrir referencia
            f.write(f"open {ref_file}\n")

            # Abrir y alinear los demás modelos a #1
            model_counter = 2  # #1 = referencia
            for output in outputs:
                file_path = os.path.abspath(getattr(self.protocol, output).getFileName())
                if file_path == ref_file:
                    continue
                f.write(f"open {file_path}\n")
                f.write(f"match #{model_counter} to #1\n")
                model_counter += 1

            # Aplicar colores **después de abrir y alinear todos**
            f.write("color byattribute occupancy palette paegreen\n")
            f.write("key darkgreen:low green: lightgreen: white:high\n")
            f.write("view orient\n")

        Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")