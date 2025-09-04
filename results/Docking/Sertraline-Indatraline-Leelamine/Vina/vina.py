import os
import subprocess

ligands = ["sertraline.pdbqt", "indatraline.pdbqt", "leelamine.pdbqt"]
receptor = "centroid1_all_atom.pdbqt"
outdir = "results"

center = {"x": 44.695, "y": 55.169, "z": 91.546}
size = {"x": 126, "y": 126, "z": 126}

exhaustiveness_list = [8, 16, 32]
replicates = 3

os.makedirs(outdir, exist_ok=True)

for lig in ligands:
    base = os.path.splitext(lig)[0]
    for exh in exhaustiveness_list:
        for r in range(1, replicates + 1):
            out_file = os.path.join(outdir, f"{base}_exh{exh}_r{r}.pdbqt")
            log_file = os.path.join(outdir, f"{base}_exh{exh}_r{r}.log")
            cmd = [
                "vina",
                "--receptor",
                receptor,
                "--ligand",
                lig,
                "--center_x",
                str(center["x"]),
                "--center_y",
                str(center["y"]),
                "--center_z",
                str(center["z"]),
                "--size_x",
                str(size["x"]),
                "--size_y",
                str(size["y"]),
                "--size_z",
                str(size["z"]),
                "--exhaustiveness",
                str(exh),
                "--out",
                out_file,
            ]
            print(f"Lancement: {lig} | Exhaustiveness={exh} | Réplicat={r}")
            subprocess.run(cmd)
