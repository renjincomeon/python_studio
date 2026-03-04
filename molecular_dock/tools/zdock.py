import os
import subprocess
import argparse

# 配置常量
DOCKER_IMAGE = "zdock302"
CONTAINER_WORKDIR = "/zdock"
WORKDIR = "/workdir"

def get_relative_to_cwd(file_path):
    """
    将绝对路径转换为相对于当前工作目录的路径。
    Docker 挂载的是 cwd 到 /zdock，所以容器内只需要相对路径。
    """
    cwd = os.getcwd()
    abs_path = os.path.abspath(file_path)
    
    if not abs_path.startswith(cwd):
        raise ValueError(f"文件 {file_path} 不在当前工作目录 {cwd} 下，Docker 无法访问。")
    
    # 获取相对路径 (例如: pdb/file.pdb)
    rel_path = os.path.relpath(abs_path, cwd)
    return rel_path

def run_pdb_mark_sur(pdb_file, config=None):
    """
    在 Docker 中运行 mark_sur
    """

    # 1️⃣ 宿主机绝对路径
    abs_input = os.path.abspath(pdb_file)
    input_dir = os.path.dirname(abs_input)
    base_name = os.path.basename(abs_input)

    # 2️⃣ 输出文件名
    if base_name.endswith(".pdb"):
        out_name = base_name[:-4] + "_m.pdb"
    else:
        out_name = base_name + "_m.pdb"

    abs_output = os.path.join(input_dir, out_name)

    # 3️⃣ 容器内部统一工作目录
    container_workdir = "/work"

    # 4️⃣ 构建 docker 命令
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{config['software']['zdock']['home']}:/zdock",  # 软件目录
        "-v", f"{input_dir}:{container_workdir}",               # 输入文件目录
        "-w", container_workdir,
        DOCKER_IMAGE,
        "/zdock/mark_sur",   # 容器内可执行文件
        base_name,           # 容器内输入文件
        out_name             # 容器内输出文件
    ]

    print(f"[mark_sur] Running: {' '.join(cmd)}")

    try:
        subprocess.run(cmd, check=True)
        print(f"[mark_sur] Success. Output: {abs_output}")
        return abs_output
    except subprocess.CalledProcessError as e:
        print(f"[mark_sur] Failed with error code {e.returncode}")
        raise

def run_zdock(receptor, ligand, outdir, config=None):
    """
    在 Docker 中运行 zdock
    """
    os.makedirs(outdir, exist_ok=True)
    
    # 1. 准备路径 (全部转为相对路径)
    rel_receptor = get_relative_to_cwd(receptor)
    rel_ligand = get_relative_to_cwd(ligand)
    
    # 输出目录的相对路径
    rel_outdir = get_relative_to_cwd(outdir)
    rel_output_file = os.path.join(rel_outdir, "zdock.out")

    # 2. 构建 Docker 命令
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{config['software']['zdock']['home']}:{CONTAINER_WORKDIR}",
         "-v", f"{dir_name}:{WORKDIR}",
        "-w", WORKDIR,
        DOCKER_IMAGE,
        f"{CONTAINER_WORKDIR}/zdock",               # 容器内执行的命令
        "-R", rel_receptor,
        "-L", rel_ligand,
        "-o", rel_output_file
    ]

    print(f"[ZDOCK] Running: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd)
        print(f"[ZDOCK] Success. Output saved to {os.path.join(outdir, 'zdock.out')}")
    except subprocess.CalledProcessError as e:
        print(f"[ZDOCK] Failed with error code {e.returncode}")
        raise
