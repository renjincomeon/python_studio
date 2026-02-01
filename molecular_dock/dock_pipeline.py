#!/usr/bin/env python3
import argparse
import yaml
import pprint
from tools.zdock import run_zdock, run_pdb_mark_sur
# from tools.hdock import run_hdock
# from tools.haddock3 import run_haddock3

def load_config(config_file="config.yaml"):
    with open(config_file, "r") as f:
        return yaml.safe_load(f)

def main():
    config = load_config("config.yaml")
    parser = argparse.ArgumentParser(
        description="Unified docking pipeline for ZDOCK / HDOCK / HADDOCK3"
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

    args = parser.parse_args()

    if "zdock" in args.methods:
        receptor_m = run_pdb_mark_sur(pdb_file=args.receptor, config=config)
        ligand_m = run_pdb_mark_sur(pdb_file=args.ligand, config=config)
        run_zdock(receptor_m, ligand_m, args.outdir + "/zdock", config)

    # if "hdock" in args.methods:
    #     run_hdock(args.receptor, args.ligand, args.outdir + "/hdock")

    # if "haddock3" in args.methods:
    #     run_haddock3(args.receptor, args.ligand, args.outdir + "/haddock3")


if __name__ == "__main__":
    # config = load_config("config.yaml")
    # pprint.pprint(config)
    main()
