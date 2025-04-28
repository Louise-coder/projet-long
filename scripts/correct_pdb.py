"""Correct a PDB file by removing H atoms and renumbering residues."""

from pathlib import Path
import os

from Bio.PDB import PDBIO, Chain, Model, PDBParser, Residue

PRED_TO_EXP = {261: 283, 609: 631, 757: 779, 1213: 1235}
WORKING_DIR = "scripts"
INPUT_FILE = f"{WORKING_DIR}/to_correct.pdb"
OUTPUT_FILE = f"{WORKING_DIR}/corrected.pdb"


def correct_pdb(pred_start: int, exp_start: int) -> None:
    """Corrects the PDB file by removing H atoms and renumbering residues.

    This function reads the PDB file line by line, removes H atoms, and
    renumbers the residues starting from the specified residue number. The
    corrected PDB file is named "corrected.pdb".

    Parameters
    ----------
    pred_start : int
        Starting residue number for the predicted PDB file.
    exp_start : int
        Starting residue number for the corrected PDB file.

    Notes
    -----
    [22:26] is the column for the residue number.
    [77:78] is the column for the atom name.

    """
    res_num = exp_start
    old_res_num = str(pred_start)
    with Path.open(INPUT_FILE) as f_to_corr, Path.open(OUTPUT_FILE, "w") as f_corr:
        for line in f_to_corr:
            if line[77:78] == "H":  # ignore H atoms
                continue
            new_res_num = str(line[22:26]).strip()
            if new_res_num != old_res_num:  # verify residue number
                res_num += 1
                old_res_num = new_res_num
            if len(new_res_num) == 3:
                f_corr.write(f"{line[:23]}{res_num}{line[26:72]}A{line[73:]}")
            else:
                f_corr.write(f"{line[:22]}{res_num}{line[26:72]}A{line[73:]}")


def correct_pdb_biopython(exp_start: int) -> None:
    """Corrects the PDB file by removing H atoms and renumbering residues.

    The function uses Biopython to parse the PDB file and create a new model
    with corrected residue numbers. The new model is saved to a file named
    "corrected.pdb".

    Parameters
    ----------
    exp_start : int
        Starting residue number for the corrected PDB file.

    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("model", INPUT_FILE)
    model = structure[0]
    new_model = Model.Model(0)
    res_num = exp_start
    new_chain = Chain.Chain("A")
    for chain in model:
        for residue in chain:
            new_residue = Residue.Residue(
                (" ", res_num, " "),
                residue.get_resname(),
                residue.get_segid(),
            )
            for atom in residue:
                if atom.element == "H":
                    continue
                new_residue.add(atom.copy())
            new_chain.add(new_residue)
            res_num += 1
        new_model.add(new_chain)
    io = PDBIO()
    io.set_structure(new_model)
    io.save(OUTPUT_FILE)


if __name__ == "__main__":
    if not os.path.exists(INPUT_FILE):
        print(f"The PDB file to correct {INPUT_FILE} does not exist.")
        exit(1)
    try:
        PRED_START = next(
            int(line[22:26]) for line in open(INPUT_FILE) if line.startswith(("ATOM"))
        )
    except StopIteration:
        print("No ATOM line found in the PDB file.")
        exit(1)
    correct_pdb(PRED_START, PRED_TO_EXP[PRED_START])
    # correct_pdb_biopython(PRED_TO_EXP[PRED_START])
