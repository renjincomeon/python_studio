#!/usr/bin/env python3
import argparse
import yaml
import pprint
import os
import sys
from tools.zdock import run_zdock, run_pdb_mark_sur
# from tools.hdock import run_hdock
# from tools.haddock3 import run_haddock3

def load_config(config_file):
    """加载配置文件，增加文件存在性检查以提供友好报错"""
    if not os.path.exists(config_file):
        print(f"Error: 配置文件 '{config_file}' 不存在。请使用 --config 指定正确的路径。")
        sys.exit(1)
    
    try:
        with open(config_file, "r", encoding='utf-8') as f:
            return yaml.safe_load(f)
    except yaml.YAMLError as e:
        print(f"Error: 解析配置文件 '{config_file}' 时出错: {e}")
        sys.exit(1)

def main():
    # 先初始化 parser
    parser = argparse.ArgumentParser(
        description="Unified docking pipeline for ZDOCK / HDOCK / HADDOCK3"
    )

    # 【新增】添加 config 文件参数
    parser.add_argument(
        "--config", 
        default="config.yaml", 
        help="Path to the configuration YAML file (default: config.yaml)"
    )

    parser.add_argument("--receptor", required=True, help="Receptor PDB file")
    parser.add_argument("--ligand", required=True, help="Ligand PDB file")
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=["zdock", "hdock", "haddock3"],
        default=["zdock", "hdock", "haddock3"],
        help="Docking methods to run"
    )
    parser.add_argument("--outdir", default="results", help="Output directory")

    # 解析参数
    args = parser.parse_args()

    # 【修改】在解析参数后，使用指定的路径加载配置
    config = load_config(args.config)

    if "zdock" in args.methods:
        # 确保输出目录存在 (可选优化)
        os.makedirs(args.outdir + "/zdock", exist_ok=True)
        
        receptor_m = run_pdb_mark_sur(pdb_file=args.receptor, config=config)
        ligand_m = run_pdb_mark_sur(pdb_file=args.ligand, config=config)
        run_zdock(receptor_m, ligand_m, args.outdir + "/zdock", config)

    # if "hdock" in args.methods:
    #     os.makedirs(args.outdir + "/hdock", exist_ok=True)
    #     run_hdock(args.receptor, args.ligand, args.outdir + "/hdock")

    # if "haddock3" in args.methods:
    #     os.makedirs(args.outdir + "/haddock3", exist_ok=True)
    #     run_haddock3(args.receptor, args.ligand, args.outdir + "/haddock3")


if __name__ == "__main__":
    # config = load_config("config.yaml") # 已移至 main 函数内部处理
    # pprint.pprint(config)
    main()
