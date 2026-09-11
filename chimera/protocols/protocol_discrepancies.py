# **************************************************************************
# *
# * Authors:   Javier Sanchez
# *             Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import subprocess
import re
import shutil
import pyworkflow.utils as pwutils
from pyworkflow.protocol import MultiPointerParam, params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import AtomicStructHandler

from pwem.objects import SetOfAtomStructs, AtomStruct
from pyworkflow.protocol import String

from chimera.protocols.protocol_base import ChimeraProtBase


class ChimeraProtDiscrepancies(EMProtocol):
    """
    Protocol to find atom discrepancies of all atomic models versus all of the rest.
    """
    _label = 'find discrepancies'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        Params:
            form: this is the form to be populated with sections and params.
        """
        form.addSection(label=Message.LABEL_INPUT)

        form.addParam('reference', params.PointerParam,
                      pointerClass='AtomStruct',
                      label="Reference structure: ",
                      help='Select the reference AtomStruct.')
        form.addParam('structures', MultiPointerParam, pointerClass="AtomStruct",
                      label='Atomic structure', important=True,
                      help='Select the set of atomic structures to be aligned and analyzed.')

    # --------------------------- INSERT steps functions ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertStep)
        self._insertFunctionStep(self.create_chimerax_script)
        self._insertFunctionStep(self.run_chimerax_script_step)
        self._insertFunctionStep(self.create_folders)
        self._insertFunctionStep(self.add_rmsd)
        self._insertFunctionStep(self.compute_mean_rmsd)
        self._insertFunctionStep(self.final_models)

    # --------------------------- STEPS functions -----------------------------
    def convertStep(self):
        self.extra_files = []
        self.model_map = {}

        # --- reference first ---
        ref_file = self.reference.get().getFileName()
        ref_origin = self.reference.get().getAttributeValue('origin')
        ref_ext = os.path.splitext(ref_file)[1]
        if ref_origin is None:
            ref_internal = f"model_00{ref_ext}"
        else:
            ref_internal = f"model_00_{ref_origin}{ref_ext}"
        ref_dest = self._getExtraPath(ref_internal)

        pwutils.createLink(ref_file, ref_dest)

        if not os.path.exists(ref_dest):
            raise Exception(f"Failed to create link for {ref_file}")

        self.extra_files.append(ref_internal)
        self.model_map[ref_internal] = os.path.basename(ref_file)

        # --- other structures ---
        for i, atomstruct in enumerate(self.structures, start=1):
            ori_file = atomstruct.get().getFileName()
            origin = atomstruct.get().getAttributeValue('origin')
            ext = os.path.splitext(ori_file)[1]
            if origin is None:
                internal_name = f"model_{i:02d}{ext}"
            else:
                internal_name = f"model_{i:02d}_{origin}{ext}"
            dest_path = self._getExtraPath(internal_name)

            pwutils.createLink(ori_file, dest_path)

            if not os.path.exists(dest_path):
                raise Exception(f"Failed to create link for {ori_file}")

            self.extra_files.append(internal_name)
            self.model_map[internal_name] = os.path.basename(ori_file)

        print(f"Conversion completed with unique names: {self.extra_files}")

    def getModelsChainsStep(self, fileName):
        """ Returns (1) list with the information
           {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
           (2) list with residues, position and chain (modelsFirstResidue)"""
        structureHandler = AtomicStructHandler()
        structureHandler.read(fileName)
        structureHandler.getStructure()
        return structureHandler.getModelsChains()

    def detect_molecule_type(self, residues):
        protein_aas = {
            "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
        }

        dna_rna = {
            "A", "C", "G", "T", "U",
            "DA", "DC", "DG", "DT", "DU", "DN", "DX"
        }

        ignored_residues = {"CL", "HOH", "NA", "MG", "K"}
        filtered_res = [r for r in residues if r not in ignored_residues]
        res_set = set(filtered_res)

        if res_set.issubset(protein_aas):
            return "protein"
        if res_set.issubset(dna_rna):
            return "nucleic"

        return "unknown"

    THREE_TO_ONE = {
        # Protein
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",

        # Nucleic acids
        "DA": "A",
        "DC": "C",
        "DG": "G",
        "DT": "T",
        "DU": "U",
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "U": "U",
    }
    def create_chimerax_script(self):
        import difflib

        project_path = self.getProject().getPath()
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        os.makedirs(output_path, exist_ok=True)
        save_path = os.path.join(project_path, output_path)

        chimerax_script = ""

        for i, model_file in enumerate(self.extra_files):
            chimerax_script += f"open {model_file}\n"

        rmsd_counter = 1
        ref_index = 1  # reference is always model 1

        ref_models, ref_residues = self.getModelsChainsStep(
            self.reference.get().getFileName()
        )
        ref_chain_dict = list(ref_models.values())[0]
        ref_chains = list(ref_chain_dict.keys())

        ref_seq_map = {}
        for ch in ref_chains:
            ref_chain_res = ref_residues[0].get(ch, [])
            ref_seq_map[ch] = [res[1] for res in ref_chain_res]

        self.chain_pairs = {}
        for i in range(len(self.extra_files)):
            model1 = os.path.splitext(self.extra_files[0])[0]
            model2 = os.path.splitext(os.path.basename(self.extra_files[i]))[0]
            if model1 == model2:
                continue

            mod_path = self._getExtraPath(self.extra_files[i])

            mod_models, mod_residues = self.getModelsChainsStep(mod_path)
            mod_chain_dict = list(mod_models.values())[0]
            mod_chains = list(mod_chain_dict.keys())

            mod_seq_map = {}
            for ch in mod_chains:
                mod_chain_res = mod_residues[0].get(ch, [])
                mod_seq_map[ch] = [res[1] for res in mod_chain_res]

            used_mod_chains = set()
            chain_pairs = []

            for ref_ch, ref_seq in ref_seq_map.items():

                best_score = 0.0
                best_mod_ch = None

                for mod_ch, mod_seq in mod_seq_map.items():

                    if mod_ch in used_mod_chains:
                        continue

                    if len(ref_seq) == 0 or len(mod_seq) == 0:
                        continue

                    from Bio import pairwise2
                    # sequence similarity (robust + simple)
                    # Convert residue names to one-letter sequences
                    three_to_one = {
                        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
                        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
                        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
                        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
                        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",

                        # DNA/RNA
                        "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
                        "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U"
                    }

                    ref_seq_1letter = ''.join(
                        three_to_one.get(res, res) for res in ref_seq
                    )

                    mod_seq_1letter = ''.join(
                        three_to_one.get(res, res) for res in mod_seq
                    )

                    # Global sequence alignment
                    alignment = pairwise2.align.globalxx(
                        ref_seq_1letter,
                        mod_seq_1letter,
                        one_alignment_only=True
                    )[0]

                    aligned_ref = alignment.seqA
                    aligned_mod = alignment.seqB

                    matches = sum(
                        a == b
                        for a, b in zip(aligned_ref, aligned_mod)
                        if a != "-" and b != "-"
                    )

                    aligned_positions = sum(
                        a != "-" and b != "-"
                        for a, b in zip(aligned_ref, aligned_mod)
                    )

                    score = matches / max(len(ref_seq_1letter), len(mod_seq_1letter))

                    print(
                        f"DEBUG SEQ MATCH: {model1} {ref_ch} ({len(ref_seq)}) "
                        f"vs {model2} {mod_ch} ({len(mod_seq)}) = {score:.4f}"
                    )

                    print(f"DEBUG 1LETTER: {ref_seq_1letter[:30]} vs {mod_seq_1letter[:30]}")

                    print(
                        f"DEBUG 1LETTER: {ref_seq_1letter[:30]} "
                        f"vs {mod_seq_1letter[:30]}"
                    )
                    if score > best_score:
                        best_score = score
                        best_mod_ch = mod_ch
                    print(f"DEBUG {model1} chain {ref_ch} FIRST 30:")
                    print(ref_seq[:30])

                    print(f"DEBUG {model2} chain {mod_ch} FIRST 30:")
                    print(mod_seq[:30])

                # accept only meaningful biological matches
                print(
                    f"DEBUG BEST MATCH: {model1} {ref_ch} -> "
                    f"{best_mod_ch} score={best_score:.4f}"
                )

                if best_mod_ch is not None and best_score >= 0.3:
                    chain_pairs.append((ref_ch, best_mod_ch))
                    used_mod_chains.add(best_mod_ch)

            if chain_pairs:
                for ch, ch2 in chain_pairs:
                    self.chain_pairs[(model1, model2, ch)] = ch2
                    ref_chain_res = ref_residues[0].get(ch, [])
                    ref_residue_names = [res[1] for res in ref_chain_res]

                    ref_type = self.detect_molecule_type(ref_residue_names)

                    if ref_type == "protein":
                        chimerax_script += (
                            f"matchmaker #{ref_index}/{ch} to #{i + 1}/{ch2} showAlignment true\n"
                        )
                        chimerax_script += "setattr a occupancy 1111.11\n"
                        chimerax_script += (
                            f"sequence header {rmsd_counter} rmsd save "
                            f"{save_path}/rmsd_{model1}_{model2}_chain_{ch}.txt\n"
                        )
                        chimerax_script += (
                            f"save {save_path}/fasta_{model1}_{model2}_chain_{ch}.fasta "
                            f"format fasta alignment {rmsd_counter}\n"
                        )
                        rmsd_counter += 1

                    elif ref_type == "nucleic":
                        chimerax_script += (
                            f"matchmaker #{ref_index}/{ch} to #{i + 1}/{ch2} showAlignment true\n"
                        )
                        chimerax_script += (
                            f"save {save_path}/fasta_{model1}_{model2}_chain_{ch}.fasta "
                            f"format fasta alignment {rmsd_counter}\n"
                        )
                        rmsd_counter += 1

            else:
                # fallback: whole-structure alignment
                chimerax_script += f"matchmaker #{ref_index} to #{i + 1} showAlignment true\n"
                chimerax_script += "setattr a occupancy 1111.11\n"
                chimerax_script += (
                    f"sequence header {rmsd_counter} rmsd save "
                    f"{save_path}/rmsd_{model1}_{model2}_chain_Whole.txt\n"
                )
                chimerax_script += (
                    f"save {save_path}/fasta_{model1}_{model2}.fasta "
                    f"format fasta alignment {rmsd_counter}\n"
                )
                rmsd_counter += 1

        for i, model_file in enumerate(self.extra_files):
            original_name = os.path.splitext(os.path.basename(model_file))[0]
            chimerax_script += f"save {save_path}/align_{original_name}.cif models #{i + 1}\n"

        chimerax_script += "exit\n"

        # save script
        script_path = os.path.join(output_path, 'chimerax_script.cxc')
        with open(script_path, 'w') as script_file:
            script_file.write(chimerax_script)

        print(f"ChimeraX script created at {script_path}")
        return script_path

    def calculate_manual_rmsd(self, cif_ref, cif_mod,
                              ref_chain_id, mod_chain_id, mol_type):
        import numpy as np

        atom = "P" if mol_type == "nucleic" else "CA"

        handler = AtomicStructHandler()

        handler.read(cif_ref)
        ref_struct = handler.getStructure()[0]

        handler.read(cif_mod)
        mod_struct = handler.getStructure()[0]
        print("DEBUG REFERENCE CHAINS:")
        for chain in ref_struct:
            print(f"  chain {chain.id}: {len(list(chain))} residues")

        print("DEBUG MODEL CHAINS:")
        for chain in mod_struct:
            print(f"  chain {chain.id}: {len(list(chain))} residues")

        print(
            f"DEBUG MANUAL RMSD: "
            f"{ref_chain_id} -> {mod_chain_id} "
            f"({mol_type})"
        )

        if ref_chain_id not in ref_struct:
            print(f"Reference chain {ref_chain_id} not found")
            return {}

        if mod_chain_id not in mod_struct:
            print(f"Model chain {mod_chain_id} not found")
            return {}

        ref_chain = list(ref_struct[ref_chain_id])
        mod_chain = list(mod_struct[mod_chain_id])
        print(f"DEBUG {cif_ref}: chain {ref_chain_id} -> {len(ref_chain)} residues")
        print(f"DEBUG {cif_mod}: chain {mod_chain_id} -> {len(mod_chain)} residues")
        print(f"DEBUG first ref residues: {[r.id for r in ref_chain[:10]]}")
        print(f"DEBUG last ref residues: {[r.id for r in ref_chain[-10:]]}")

        ref_coords = []
        mod_coords = []
        residue_ids = []

        n = min(len(ref_chain), len(mod_chain))

        for i in range(n):
            r1 = ref_chain[i]
            r2 = mod_chain[i]

            if atom in r1 and atom in r2:
                ref_coords.append(r1[atom].get_coord())
                mod_coords.append(r2[atom].get_coord())
                residue_ids.append(r1.id[1])

        if len(ref_coords) < 3:
            return {}

        P = np.array(ref_coords)
        Q = np.array(mod_coords)

        Pc = P - P.mean(axis=0)
        Qc = Q - Q.mean(axis=0)

        C = Pc.T @ Qc
        V, S, Wt = np.linalg.svd(C)

        if np.linalg.det(V @ Wt) < 0:
            V[:, -1] *= -1

        U = V @ Wt

        P_aligned = Pc @ U

        rmsd_per_res = {}

        for res_id, p, q in zip(residue_ids, P_aligned, Qc):
            rmsd_per_res[res_id] = float(np.linalg.norm(p - q))

        return rmsd_per_res

    def run_chimerax_script(self, script_path, output_log_path):
        if not os.path.exists(script_path):
            raise Exception(f"ChimeraX script not found at {script_path}")

        # Get the path for the Chimera executable - Version 1.11.1
        self.scipion_path = os.environ.get('SCIPION_HOME', None)
        if self.scipion_path is None:
            raise Exception(
                "SCIPION_HOME environment variable is not set. Please set it to the Scipion installation directory.")
        #self.chimerax_executable = f"flatpak run edu.ucsf.rbvi.ChimeraX"
        self.chimerax_executable = os.path.join(self.scipion_path, 'software/em/chimerax-1.6.1/bin/ChimeraX')

        # Run the ChimeraX script and capture the output
        print(f"Chimera: {self.chimerax_executable}")
        #result = subprocess.run(f"{self.chimerax_executable} --nogui {script_path}",
        #                        shell=True, capture_output=True, text=True)
        result = subprocess.run([self.chimerax_executable, '--nogui', script_path], capture_output=True, text=True)

        # Save the output log to a file
        with open(output_log_path, 'w') as log_file:
            log_file.write(result.stdout)

        if result.returncode != 0:
            with open(output_log_path, 'r') as log_file:
                log_contents = log_file.read()
            raise Exception(
                f"ChimeraX script failed with return code {result.returncode}. See log for details:\n{log_contents}")

        print(f"ChimeraX script executed successfully. Log saved at {output_log_path}")

    def run_chimerax_script_step(self):
        script_path = os.path.join(self.getWorkingDir(), 'extra', 'chimerax_script.cxc')
        output_log_path = os.path.join(self.getWorkingDir(), 'extra', 'chimerax_output.log')

        self.run_chimerax_script(script_path, output_log_path)
        print(f"ChimeraX script executed. Log saved at {output_log_path}")

    def create_folders(self):
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        for fasta_file in os.listdir(output_path):
            if fasta_file.startswith('fasta_') and fasta_file.endswith('.fasta'):
                pair = fasta_file.replace("fasta_", "").replace(".fasta", "")
                pair_no_chain = pair.split("_chain_")[0]
                model1, model2 = pair_no_chain.split("_model_")
                model2 = "model_" + model2

                folder_name = f"{model1}_{model2}"
                folder_path = os.path.join(output_path, folder_name)
                os.makedirs(folder_path, exist_ok=True)

                # Link the FASTA / RMSD / CIF files
                src_fasta_file = os.path.join(output_path, fasta_file)
                dest_fasta_file = os.path.join(folder_path, fasta_file)
                pwutils.createLink(src_fasta_file, dest_fasta_file)

                rmsd_files = [
                    f for f in os.listdir(output_path)
                    if f.startswith(f"rmsd_{model1}_{model2}_chain_")
                ]
                for rmsd_file in rmsd_files:
                    src_rmsd_file = os.path.join(output_path, rmsd_file)
                    dest_rmsd_file = os.path.join(folder_path, rmsd_file)
                    pwutils.createLink(src_rmsd_file, dest_rmsd_file)

                cif_file = f"align_{model1}.cif"
                src_cif_file = os.path.join(output_path, cif_file)
                dest_cif_file = os.path.join(folder_path, cif_file)
                pwutils.createLink(src_cif_file, dest_cif_file)
                out_model1_file = os.path.join(folder_path, f"out_{model1}.cif")
                pwutils.copyFile(src_cif_file, out_model1_file)
                cif_file = f"align_{model2}.cif"
                src_cif_file = os.path.join(output_path, cif_file)
                dest_cif_file = os.path.join(folder_path, cif_file)
                pwutils.createLink(src_cif_file, dest_cif_file)
                out_model2_file = os.path.join(folder_path, f"out_{model2}.cif")
                pwutils.copyFile(src_cif_file, out_model2_file)

        for folder in os.listdir(output_path):
            folder_path = os.path.join(output_path, folder)
            if os.path.isdir(folder_path) and "_" in folder:
                model1, model2 = folder.rsplit("_model_", 1)
                model2 = "model_" + model2
                out_model1_file = os.path.join(folder_path, f"out_{model1}.cif")
                out_model2_file = os.path.join(folder_path, f"out_{model2}.cif")

                text_to_append = "\nloop_\n_scipion_attributes.name\n_scipion_attributes.recipient\n_scipion_attributes.specifier\n_scipion_attributes.value\n"

                with open(out_model1_file, 'a') as file1:
                    file1.write(text_to_append)

                with open(out_model2_file, 'a') as file2:
                    file2.write(text_to_append)

    def add_rmsd(self):
        output_path = os.path.join(self.getWorkingDir(), 'extra')
        self.aa_rmsd = []
        self.occ_position = [[], []]

        for folder in os.listdir(output_path):
            folder_path = os.path.join(output_path, folder)

            if os.path.isdir(folder_path) and "_" in folder:
                model1, model2 = folder.rsplit("_model_", 1)
                model2 = "model_" + model2

                fasta_files = [
                    f for f in os.listdir(folder_path)
                    if f.startswith(f"fasta_{model1}_{model2}_chain_")
                ]

                for fasta_file in fasta_files:

                    chain = fasta_file.split("_chain_")[-1].replace(".fasta", "")

                    fasta_path = os.path.join(folder_path, fasta_file)
                    rmsd_path = os.path.join(
                        folder_path,
                        f"rmsd_{model1}_{model2}_chain_{chain}.txt"
                    )

                    rmsd_valid = False

                    if os.path.exists(rmsd_path):
                        with open(rmsd_path, 'r') as f:
                            rmsd_lines = [line.strip() for line in f if line.strip()]

                        # Header + at least one RMSD value
                        rmsd_valid = len(rmsd_lines) > 1

                    if not rmsd_valid:
                        cif_ref = os.path.join(folder_path, f"align_{model1}.cif")
                        cif_mod = os.path.join(folder_path, f"align_{model2}.cif")

                        # Get the chain corresponding to the reference chain
                        mod_chain = self.chain_pairs.get(
                            (model1, model2, chain)
                        )

                        if mod_chain is None:
                            print(
                                f"WARNING: No chain mapping found for "
                                f"{model1} chain {chain} -> {model2}"
                            )
                            manual_rmsd = {}
                        else:
                            # Determine whether this chain is protein or nucleic acid
                            ref_residues = self.getModelsChainsStep(cif_ref)[1][0]
                            ref_chain_res = ref_residues.get(chain, [])
                            ref_residue_names = [res[1] for res in ref_chain_res]

                            mol_type = self.detect_molecule_type(ref_residue_names)

                            print(
                                f"DEBUG FALLBACK: {model1} chain {chain} "
                                f"-> {model2} chain {mod_chain} | type={mol_type}"
                            )

                            manual_rmsd = self.calculate_manual_rmsd(
                                cif_ref,
                                cif_mod,
                                chain,
                                mod_chain,
                                mol_type
                            )

                        with open(rmsd_path, 'w') as f:
                            f.write("residue:rmsd\n")
                            for res, val in manual_rmsd.items():
                                f.write(f"{res}:{val}\n")

                        root_extra = self._getExtraPath()
                        root_rmsd_path = os.path.join(
                            root_extra,
                            os.path.basename(rmsd_path)
                        )

                        if not os.path.exists(root_rmsd_path):
                            pwutils.createLink(rmsd_path, root_rmsd_path)

                    # ---------------- FASTA ----------------
                    mat1, mat2 = [], []
                    with open(fasta_path, 'r') as f:
                        lines = f.readlines()
                        reading_model1, reading_model2 = False, False

                        for line in lines:
                            if line.startswith(f">{model1}"):
                                reading_model1 = True
                                reading_model2 = False
                                continue
                            elif line.startswith(f">{model2}"):
                                reading_model1 = False
                                reading_model2 = True
                                continue
                            elif line.startswith(">"):
                                reading_model1 = reading_model2 = False

                            if reading_model1:
                                mat1.extend([0 if c == '.' else 1 for c in line.strip()])
                            elif reading_model2:
                                mat2.extend([0 if c == '.' else 1 for c in line.strip()])

                    mat_rmsd = []
                    with open(rmsd_path, 'r') as f:
                        lines = f.readlines()[1:]

                        for line in lines:
                            parts = line.split(":")
                            if len(parts) != 2:
                                continue

                            try:
                                Aa = int(parts[0].strip())
                                rmsd_value = float(parts[1].strip()) if parts[1].strip() != "None" else 250.0

                                while len(mat_rmsd) <= Aa:
                                    mat_rmsd.append(0)

                                mat_rmsd[Aa] += rmsd_value
                            except ValueError:
                                continue

                    def build_positions(mat):
                        pos_align, pos_res = [], []
                        pos = 0
                        for i, v in enumerate(mat):
                            if v != 0:
                                pos_align.append(i)
                                pos += 1
                                pos_res.append(pos)
                        return pos_align, pos_res

                    mat_pos1_align, mat_pos1_resnum = build_positions(mat1)
                    mat_pos2_align, mat_pos2_resnum = build_positions(mat2)

                    def build_occ(pos_align, pos_res):
                        occ_dict = {}
                        aa_occ = [[], []]

                        for i, resnum in enumerate(pos_res):
                            align_idx = pos_align[i]
                            rmsd_value = mat_rmsd[align_idx] if align_idx < len(mat_rmsd) else 250.0

                            occ_dict[resnum] = rmsd_value
                            aa_occ[0].append(resnum)
                            aa_occ[1].append(rmsd_value)

                        return occ_dict, aa_occ

                    occupancy1_dict, aa_occ1 = build_occ(mat_pos1_align, mat_pos1_resnum)
                    occupancy2_dict, aa_occ2 = build_occ(mat_pos2_align, mat_pos2_resnum)

                    self.aa_rmsd.append({f"{model1}_chain_{chain}": aa_occ1})
                    self.aa_rmsd.append({f"{model2}_chain_{chain}": aa_occ2})

                    out_model1_file = os.path.join(folder_path, f"out_{model1}.cif")
                    out_model2_file = os.path.join(folder_path, f"out_{model2}.cif")

                    with open(out_model1_file, 'r') as f1, open(out_model2_file, 'r') as f2:
                        lines1 = f1.readlines()
                        lines2 = f2.readlines()

                    subs = []

                    def process(lines, occ_dict, out):
                        for line in lines:
                            if line.startswith('ATOM'):
                                clean = re.sub(r'\s+', ' ', line).strip()
                                cols = clean.split()

                                if '1111.11' in cols:
                                    subs.append(cols.index('1111.11'))

                                resnum = int(cols[8]) if cols[8].isdigit() else None

                                if resnum:
                                    val = round(occ_dict.get(resnum, 0.0), 2)
                                    line = re.sub(r'1111\.11', str(val), line)

                            out.write(line)

                    with open(out_model1_file, 'w') as o1, open(out_model2_file, 'w') as o2:
                        process(lines1, occupancy1_dict, o1)
                        process(lines2, occupancy2_dict, o2)

                    subs = list(set(subs))
                    self.occ_position[0] += [model1, model2]

                    if subs:
                        self.occ_position[1] += [subs[0], subs[0]]

    def compute_mean_rmsd(self):
        output_path = os.path.join(self.getWorkingDir(), 'extra')

        rmsd_by_chain = {}

        rmsd_files = [
            f for f in os.listdir(output_path)
            if f.startswith('rmsd_') and '_chain_' in f
        ]

        for rmsd_file in rmsd_files:
            chain = rmsd_file.split("_chain_")[-1].replace(".txt", "")

            file_path = os.path.join(output_path, rmsd_file)

            vals = []
            with open(file_path, 'r') as f:
                lines = f.readlines()[1:]

                for line in lines:
                    parts = line.strip().split(":")
                    if len(parts) == 2:
                        try:
                            val = float(parts[1].strip()) if parts[1].strip() != "None" else 0.0
                            vals.append(val)
                        except:
                            continue

            if chain not in rmsd_by_chain:
                rmsd_by_chain[chain] = []

            rmsd_by_chain[chain].append(vals)

        for chain, matrices in rmsd_by_chain.items():
            mean_rmsd = []

            max_len = max(len(m) for m in matrices)

            for i in range(max_len):
                vals = [m[i] for m in matrices if i < len(m)]
                mean_rmsd.append(sum(vals) / len(vals) if vals else 0.0)

            out_file = os.path.join(output_path, f"mean_rmsd_chain_{chain}.txt")

            with open(out_file, 'w') as f:
                for i, val in enumerate(mean_rmsd, start=1):
                    f.write(f"{i}: {val:.6f}\n")

            print(f"Mean RMSD chain {chain}: {out_file}")

    def final_models(self):
        import os, shutil, re
        from pwem.objects import AtomStruct

        output_path = os.path.join(self.getWorkingDir(), 'extra')
        final_output_path = os.path.join(output_path, 'FINAL-OUTPUTS')
        os.makedirs(final_output_path, exist_ok=True)

        for folder in os.listdir(output_path):
            folder_path = os.path.join(output_path, folder)

            if os.path.isdir(folder_path) and folder != 'FINAL-OUTPUTS':
                for file_name in os.listdir(folder_path):
                    if file_name.startswith('out_') and file_name.endswith('.cif'):
                        shutil.copy(
                            os.path.join(folder_path, file_name),
                            os.path.join(final_output_path, file_name)
                        )

        grouped = {}

        for entry in self.aa_rmsd:
            for key, matrix in entry.items():
                if key not in grouped:
                    grouped[key] = []
                grouped[key].append(matrix)

        averaged = {}

        for key, matrices in grouped.items():
            max_len = max(len(m[1]) for m in matrices)

            avg = []
            for i in range(max_len):
                vals = [m[1][i] for m in matrices if i < len(m[1])]
                avg.append(sum(vals) / len(vals) if vals else 0.0)

            averaged[key] = avg

        for file_name in os.listdir(final_output_path):
            if not file_name.startswith("out_") or not file_name.endswith(".cif"):
                continue

            model = file_name[4:-4]
            file_path = os.path.join(final_output_path, file_name)

            with open(file_path, 'r') as f:
                lines = f.readlines()

            updated_lines = []
            counters = {}

            for line in lines:
                if line.startswith('ATOM'):
                    clean = re.sub(r'\s+', ' ', line).strip()
                    cols = clean.split()

                    chain = cols[4] if len(cols) > 4 else None
                    resnum = int(cols[8]) if len(cols) > 8 and cols[8].isdigit() else None

                    if chain is not None and resnum is not None:
                        key = f"{model}_chain_{chain}"

                        if key in averaged:
                            if key not in counters:
                                counters[key] = 0

                            idx = counters[key]

                            if idx < len(averaged[key]):
                                val = round(averaged[key][idx], 2)
                                line = re.sub(r'1111\.11', str(val), line)

                            counters[key] += 1

                updated_lines.append(line)

            with open(file_path, 'w') as f:
                f.writelines(updated_lines)

        ref_base = os.path.splitext(self.extra_files[0])[0]

        for file_name in os.listdir(final_output_path):
            file_path = os.path.join(final_output_path, file_name)

            if os.path.isfile(file_path) and file_name.startswith("out_"):
                output = AtomStruct(filename=file_path)

                output_name = os.path.splitext(file_name)[0]
                file_base = output_name.replace("out_", "").lower()

                if file_base == ref_base:
                    self._defineOutputs(**{f'ref_{output_name}': output})
                else:
                    self._defineOutputs(**{output_name: output})

        print("Final models with per-chain RMSD generated")

