import os
import pyworkflow.viewer as pwviewer
from pwem.viewers.viewer_chimera import Chimera
from ..protocols.protocol_discrepancies import ChimeraProtDiscrepancies
import pyworkflow.protocol.params as params
import matplotlib.pyplot as plt
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
        group.addParam('showFiles', params.EnumParam, choices=['all']+self.getFileNames(),
                       label='Plot:', default=0,
                       help='Plot all files or specific one')
        group.addParam(
                        'selectedChain', params.EnumParam, choices=['all'] + self._getAvailableChains(),
                        default=0, label='Chain:',
                        help='Select chain to display')
        group.addParam('filter',
                       params.BooleanParam, default=False,
                       label='Filter residues:',
                       help='Generate and display a graph of only selected residues')
        group.addParam('residuesMin', params.IntParam,
                       label='From:', default=0, condition='filter',
                       help='From x residue')
        group.addParam('residuesMax', params.IntParam,
                       label='To:', default=1000, condition='filter',
                       help='To x residue')

    def getFileNames(self):
        outputs = self.protocol._outputs
        cleanList = []
        for out in outputs:
            clean = self._normalizeOutputName(out)
            if '00' not in clean:
                cleanList.append(clean)
        return cleanList

    def _getAvailableChains(self):
        extra_path = self.protocol._getExtraPath()
        files = [f for f in os.listdir(extra_path) if f.startswith("rmsd_") and f.endswith(".txt")]
        chains = set()
        for f in files:
            # Only count actual chains, ignore the ones labeled _Whole
            if "_chain_" in f and "_Whole" not in f:
                chain = f.split("_chain_")[-1].replace(".txt", "")
                chains.add(chain)
        return sorted(list(chains))

    def _normalizeOutputName(self, output_name):
        name = output_name.replace("out_", "").replace("ref_out_", "")
        return name

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

            if not outputs:
                return []

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
            if f.startswith("rmsd_") and f.endswith(".txt")
        ])

        choices = ['all'] + self._getAvailableChains()
        selected_chain = choices[self.selectedChain.get()]

        selectedIdx = self.showFiles.get()
        choices = ['all'] + self.getFileNames()
        selectedValue = choices[selectedIdx]
        if selectedValue != 'all':
            normalized = self._normalizeOutputName(selectedValue)
            files = [f for f in files if normalized in f]

        if not files:
            return []

        chains = {}
        for f in files:
            if "_chain_" in f:
                chain = f.split("_chain_")[-1].replace(".txt", "")
            else:
                chain = "WholeModel"
            if selected_chain != 'all' and chain != selected_chain:
                continue

            chains.setdefault(chain, []).append(f)

        def get_model_name(fname):
            base = fname.replace("rmsd_", "").replace(".txt", "")
            if "_chain_" in base:
                return base.split("_chain_")[0]
            return base


        models = sorted({get_model_name(f) for f in files})
        colors = cm.tab10.colors
        color_map = {m: colors[i % len(colors)] for i, m in enumerate(models)}

        fig, ax = plt.subplots(figsize=(10, 6))
        fig.canvas.manager.set_window_title(
            f"RMSD_Protocol_{self.protocol.getObjId()}"
        )

        current_x_offset = 0
        all_positions = []
        all_labels = []

        is_filtered = self.filter.get()
        min_r = self.residuesMin.get() if is_filtered else -float('inf')
        max_r = self.residuesMax.get() if is_filtered else float('inf')

        from matplotlib.lines import Line2D

        for chain, chain_files in sorted(chains.items()):
            chain_max_x = 0
            chain_residues_plotted = set()

            for f in chain_files:
                path = os.path.join(extra_path, f)
                model = get_model_name(f)
                color = color_map[model]

                x_vals, y_vals = [], []

                with open(path, 'r') as fh:
                    for line in fh:
                        if ':' not in line:
                            continue

                        try:
                            res_id = int(line.split(':')[0].strip())
                            rmsd_val = float(line.split(':')[1].strip())

                            if is_filtered and (res_id < min_r or res_id > max_r):
                                continue

                            plot_x = (res_id - min_r) if is_filtered else res_id

                            x_coord = plot_x + current_x_offset

                            x_vals.append(x_coord)
                            y_vals.append(rmsd_val)

                            chain_residues_plotted.add((x_coord, res_id))

                        except Exception:
                            continue

                if x_vals:
                    ax.plot(
                        x_vals,
                        y_vals,
                        color=color,
                        linewidth=1.2,
                        alpha=0.9
                    )

                    chain_max_x = max(
                        chain_max_x,
                        max(x_vals) - current_x_offset
                    )

            if not chain_residues_plotted:
                continue

            ax.axvline(
                current_x_offset - 1,
                linestyle='--',
                color='gray',
                alpha=0.2
            )

            chain_ticks = sorted(list(chain_residues_plotted))

            all_positions.extend([t[0] for t in chain_ticks])
            all_labels.extend([t[1] for t in chain_ticks])

            chain_center = current_x_offset + chain_max_x / 2

            y_top = ax.get_ylim()[1]
            y_offset = (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.03  # small padding

            ax.text(
                chain_center,
                y_top - y_offset,
                f"Chain {chain}",
                ha='center',
                va='bottom',
                fontsize=10,
                fontweight='bold',
                color='black',
                alpha=0.9
            )

            current_x_offset += chain_max_x + 5

        handles = [
            Line2D(
                [0],
                [0],
                color=color_map[m],
                lw=2,
                label=m
            )
            for m in models
        ]

        if all_positions:
            step = max(1, len(all_labels) // 30)

            ax.set_xticks(all_positions[::step])

            ax.set_xticklabels(
                all_labels[::step],
                rotation=90,
                fontsize=6
            )

        ax.set_xlabel("Residue ID")
        ax.set_ylabel("RMSD")
        ax.set_title("Filtered RMSD per residue (Grouped by Chain)")


        ax.grid(True, alpha=0.3)

        # Excellent agreement
        ax.axhline(
            y=2.0,
            color='green',
            linestyle='--',
            linewidth=0.5,
            alpha=0.8
        )
        # Good agreement
        ax.axhline(
            y=4.0,
            color='red',
            linestyle='--',
            linewidth=0.5,
            alpha=0.8
        )

        handles.extend([
            Line2D([0], [0], color='green', linestyle='--', lw=1.2,
                   label='Excellent (<2 Å)'),
            Line2D([0], [0], color='red', linestyle='--', lw=1.5,
                   label='Good (<4 Å)')
        ])
        ax.legend(handles=handles, fontsize=8)

        fig.tight_layout()

        plt.show(block=False)

        return []