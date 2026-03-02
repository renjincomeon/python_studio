import os
import subprocess

HOST_ROOT_DIR = os.getcwd()
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
    result_file = pdb_file.replace('.pdb', '_m.pdb')
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{HOST_ROOT_DIR}:/zdock -w /zdock zdock302",
        "-w", "/zdock",
        "zdock302",
        "./mark_sur",
        pdb_file,
        result_file
    ]

    print("[mark_sur] Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return result_file
