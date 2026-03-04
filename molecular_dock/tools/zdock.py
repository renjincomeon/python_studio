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
    # 1. 准备路径
    # 输入文件的相对路径 (容器内可见)
    rel_input = get_relative_to_cwd(pdb_file)
    
    # 输出文件逻辑：同目录下，文件名加 _m.pdb
    dir_name = os.path.dirname(rel_input)
    base_name = os.path.basename(rel_input)
    if base_name.endswith('.pdb'):
        out_name = base_name[:-4] + "_m.pdb"
    else:
        out_name = base_name + "_m.pdb"
    
    rel_output = os.path.join(dir_name, out_name)
    
    # 宿主机的完整输出路径 (用于返回给 Python)
    abs_output = os.path.join(os.getcwd(), rel_output)

    # 2. 构建 Docker 命令 (关键修复点)
    # 格式: docker run --rm -v <cwd>:/zdock -w /zdock <image> <command> <args>
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{config['software']['zdock']['home']}:{CONTAINER_WORKDIR}",
        "-v", f"{dir_name}:{WORKDIR}",
        "-w", WORKDIR,
        DOCKER_IMAGE,
        f"{CONTAINER_WORKDIR}/mark_sur",          # 容器内执行的命令
        rel_input,             # 容器内输入路径 (相对路径)
        rel_output             # 容器内输出路径 (相对路径)
    ]

    print(f"[mark_sur] Running: {' '.join(cmd)}")
    
    try:
        subprocess.run(cmd)
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
