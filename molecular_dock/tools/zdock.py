import os
import subprocess


def run_zdock(receptor, ligand, outdir, config):
    os.makedirs(outdir, exist_ok=True)

    cmd = [
        "zdock",
        "-R", receptor,
        "-L", ligand,
        "-o", os.path.join(outdir, "zdock.out")
    ]

    print("[ZDOCK] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

def run_pdb_mark_sur(pdb_file, config):
    soft_path = config['software']['zdock']['home']
    soft = os.path.join(soft_path, 'mark_sur')
    result_file = pdb_file.replace('.pdb', '_m.pdb')
    cmd = [
        soft,
        pdb_file,
        result_file
    ]

    print("[mark_sur] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return result_file
