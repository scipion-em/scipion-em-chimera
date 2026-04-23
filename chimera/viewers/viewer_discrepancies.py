import os
import pyworkflow.viewer as pwviewer
from pwem.viewers.viewer_chimera import Chimera
from ..protocols.protocol_discrepancies import ChimeraProtDiscrepancies
import pyworkflow.protocol.params as params


class ChimeraProtDiscrepanciesViewer(pwviewer.ProtocolViewer):
    """ Viewer for ChimeraProtDiscrepancies protocol output. """
    _label = 'viewer discrepancies'
    _targets = [ChimeraProtDiscrepancies]

    def __init__(self, **args):
        super().__init__(**args)

    def _defineParams(self, form):
        form.addSection(label='Visualization of discrepancies')
        group = form.addGroup('Open Chimera ')
        group.addParam('displayStructs',
                       params.LabelParam,
                       label='Open structures in Chimera: ',
                       help='Display the structures in Chimera GUI.')
        group = form.addGroup('Graph')
        group.addParam('displayRMSD',
                       params.LabelParam,
                       label='Graph RMSD:',
                       help='Generate and display a graph')

    def _getVisualizeDict(self):
        visDic = super()._getVisualizeDict()

        visDic.update({
            'displayStructs': self.viewChimera,
            'displayRMSD': self.viewGraphRMSD,
        })
        return visDic

    def viewChimera(self,  paramName=None):
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

        return [Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")]

    def viewGraphRMSD(self, paramName=None):
        import matplotlib.pyplot as plt
        import os

        extra_path = self.protocol._getExtraPath()
        files = sorted([
            f for f in os.listdir(extra_path)
            if f.startswith("rmsd_model_") and f.endswith(".txt")
        ])

        if not files:
            return []

        plt.figure(figsize=(10, 6))

        for f in files:
            x, y = [], []
            path = os.path.join(extra_path, f)

            with open(path, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if not line or ':' not in line:
                        continue
                    try:
                        parts = line.split(':')
                        # Only add if both sides of the colon exist
                        res_id = int(parts[0].strip())
                        rmsd_val = float(parts[1].strip())
                        x.append(res_id)
                        y.append(rmsd_val)
                    except (ValueError, IndexError):
                        continue

            # --- THE CRITICAL CHECK ---
            if len(x) != len(y):
                print(f"WARNING: Skipping {f} due to mismatched data: X={len(x)}, Y={len(y)}")
                continue

            if x and y:
                label = f.replace("rmsd_model_", "").replace(".txt", "")
                plt.plot(x, y, label=label)

        plt.xlabel("Residue")
        plt.ylabel("Mean RMSD")
        plt.title("Mean RMSD per residue")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.show()
        plt.close()

        return []