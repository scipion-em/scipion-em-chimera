import os
import pyworkflow.viewer as pwviewer
from pwem.viewers.viewer_chimera import Chimera
from ..protocols.protocol_discrepancies import ChimeraProtDiscrepancies
import pyworkflow.protocol.params as params
import matplotlib.pyplot as plt
import re
import matplotlib.cm as cm


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

    def viewChimera(self, paramName=None):
        import os

        fnCmd = self.protocol._getExtraPath("discrepancies_viewer_with_files.cxc")

        with open(fnCmd, 'w') as f:

            outputs = self.protocol._outputs
            print(outputs)

            if not outputs:
                return []

            # --------------------------
            # REFERENCE
            # --------------------------
            ref_file = None

            for output in outputs:
                obj = getattr(self.protocol, output, None)
                if obj is None:
                    continue

                file_path = os.path.abspath(obj.getFileName())

                if os.path.basename(file_path).startswith("ref_"):
                    ref_file = file_path
                    break

            if not ref_file:
                obj0 = getattr(self.protocol, outputs[0], None)
                if obj0 is None:
                    return []
                ref_file = os.path.abspath(obj0.getFileName())

            f.write(f"open {ref_file}\n")

            model_counter = 2

            for output in outputs:
                obj = getattr(self.protocol, output, None)
                if obj is None:
                    continue

                file_path = os.path.abspath(obj.getFileName())

                if file_path == ref_file:
                    continue

                f.write(f"open {file_path}\n")
                f.write(f"match #{model_counter} to #1\n")

                model_counter += 1

            f.write("color byattribute occupancy palette paegreen\n")
            f.write("key darkgreen:low green: lightgreen: white:high\n")
            f.write("view orient\n")

        return [Chimera.runProgram(Chimera.getProgram(), fnCmd + "&")]

    def viewGraphRMSD(self, paramName=None):
        extra_path = self.protocol._getExtraPath()

        files = sorted([
            f for f in os.listdir(extra_path)
            if f.startswith("rmsd_") and "_chain_" in f and f.endswith(".txt")
        ])

        if not files:
            return []

        chains = {}
        for f in files:
            chain = f.split("_chain_")[-1].replace(".txt", "")
            chains.setdefault(chain, []).append(f)

        def get_model_name(fname):
            base = fname.replace("rmsd_", "").replace(".txt", "")
            return base.split("_chain_")[0]  # TODO lo importante

        models = sorted({get_model_name(f) for f in files})

        colors = cm.tab10.colors
        color_map = {m: colors[i % len(colors)] for i, m in enumerate(models)}

        plt.figure(figsize=(10, 6))

        offset = 0
        xticks = []
        xtick_labels = []

        for chain, chain_files in sorted(chains.items()):

            max_len_chain = 0

            for f in chain_files:
                path = os.path.join(extra_path, f)

                model = get_model_name(f)
                color = color_map[model]

                x, y = [], []

                with open(path, 'r') as fh:
                    for line in fh:
                        if ':' not in line:
                            continue
                        try:
                            res_id = int(line.split(':')[0].strip())
                            rmsd_val = float(line.split(':')[1].strip())

                            x.append(res_id + offset)
                            y.append(rmsd_val)
                        except:
                            continue

                if not x:
                    continue

                plt.plot(
                    x,
                    y,
                    color=color,
                    linewidth=1.2,
                    alpha=0.9
                )

                max_len_chain = max(max_len_chain, max(x) - offset)

            plt.axvline(offset, linestyle='--', alpha=0.2)

            xticks.append(offset + max_len_chain // 2)
            xtick_labels.append(chain)

            offset += max_len_chain + 2

        handles = [
            plt.Line2D([0], [0], color=color_map[m], lw=2, label=m)
            for m in models
        ]

        all_positions = []
        all_labels = []

        offset = 0

        for chain, chain_files in sorted(chains.items()):

            max_len_chain = 0

            sample_file = chain_files[0]
            path = os.path.join(extra_path, sample_file)

            residues = []

            with open(path, 'r') as fh:
                for line in fh:
                    if ':' not in line:
                        continue
                    try:
                        res_id = int(line.split(':')[0].strip())
                        residues.append(res_id)
                    except:
                        continue

            if residues:
                positions = [r + offset for r in residues]

                all_positions.extend(positions)
                all_labels.extend(residues)

                max_len_chain = max(residues)

            offset += max_len_chain + 2

        step = max(1, len(all_labels) // 30)

        plt.xticks(all_positions[::step], all_labels[::step], rotation=90, fontsize=6)

        plt.xlabel("Residue (chains concatenated)")
        plt.ylabel("RMSD")
        plt.title("RMSD per residue grouped by chain (same color per model)")
        plt.legend(handles=handles, fontsize=8)
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()
        plt.close()

        return []