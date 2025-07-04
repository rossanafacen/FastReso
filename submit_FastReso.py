"""
Module to submit FastReso jobs 
"""
import os
from pathlib import Path
import logging
import argparse
import sys
import yaml
import subprocess
import json



def write_all_T(config, T_type):
    dT = config['physics']['dT']
    n_steps = int((config['physics'][T_type][1] - config['physics'][T_type][0])/dT) + 1
    config['physics'][T_type] = [round(config['physics'][T_type][0] + i * dT,3) for i in range(n_steps)]
    return 0


def submit_to_slurm(out_path, job, config):
    scriptname = f"slurm_{config['slurm']['name']}.sh"
    if config['physics']['do_pce'] == True:
        n_jobs = int(len(config['physics']['Tchem']))
    else:
        n_jobs = int(len(config['physics']['Tkin']))
        
    with open(Path(out_path, scriptname), "w", encoding="utf-8") as fbash:
        fbash.write("#!/bin/bash\n")
        fbash.write(f"#SBATCH --array=0-{n_jobs-1} \n")
        fbash.write(f"#SBATCH --job-name={config['slurm']['name']}\n")
        fbash.write(f"#SBATCH --time={config['slurm']['max_time']}\n")
        fbash.write(f"#SBATCH --partition={config['slurm']['partition']}\n")
        fbash.write(f"#SBATCH --chdir={config['slurm']['work_directory']}\n")
        fbash.write(f"#SBATCH --mem={config['slurm']['memory']}\n")
        fbash.write(f"#SBATCH --cpus-per-task={config['slurm']['cpus_per_job']}\n")
        fbash.write("#SBATCH -o logs/%x_%A_%a.out\n")
        fbash.write("#SBATCH -e logs/%x_%A_%a.err\n")
        fbash.write("\n")
        fbash.write(f"{job}")
        os.chmod(fbash.name, 0o755)
    
    script_path = Path(out_path, scriptname)
    subprocess.run(["sbatch", script_path], capture_output=True, cwd=out_path, text=True)
    logging.info("submitted %d jobs ", n_jobs)
    return 0


def main() -> int:

    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser(description="submit slurm jobs.")

    parser.add_argument("config", type=str, help="configuration file")
    
    args = parser.parse_args()
    config_path = Path(args.config)
    with open(config_path, "r", encoding="utf-8") as f:
        config = yaml.safe_load(f)
   
    out_path = Path(config["out_directory"])
    if out_path.is_dir():
        logging.fatal("Directory already exists. Output can not be overwritten. Please fix!")
        return 1
    
    out_path.mkdir(exist_ok=True, parents=True)

    Path(out_path, "logs").mkdir(exist_ok=True, parents=True)
    Path(out_path, "yaml").mkdir(exist_ok=True, parents=True)
    Path(out_path, "Fj").mkdir(exist_ok=True, parents=True)
    Path(out_path, "Kj").mkdir(exist_ok=True, parents=True)

    Fj_dir = Path(out_path, "Fj")
    Kj_dir = Path(out_path, "Kj")

    
    
    if config['physics']['do_pce'] == True:
        write_all_T(config, 'Tchem')
        run_Fj = config['run_pce']
        decay_list = f"{config['reversed_decay_list']} {config['decay_list']}"
        run_Kj = config['run_Kj_pce']
        
    else:
        write_all_T(config, 'Tkin')
        run_Fj = config['run_Fj']
        decay_list = config['decay_list']
        run_Kj = config['run_Kj']


    with open(Path(out_path, "configuration.json"), "w", encoding="utf-8") as fconf:
        json.dump(config, fconf, indent=4)
    
    configuration_path = Path(out_path, "configuration.json")

    job = (
        f"cp {config_path.absolute()} {out_path.absolute()}/yaml/analysis_$SLURM_JOB_ID.yaml \n"
        f"{run_Fj} {decay_list} {Fj_dir}/ {configuration_path} $SLURM_ARRAY_TASK_ID \n"
        f"{run_Kj} {Fj_dir}/ {Kj_dir}/ {configuration_path} $SLURM_ARRAY_TASK_ID \n"
    )
    
    submit_to_slurm(out_path, job, config)

    return 0


if __name__ == "__main__":
    sys.exit(main())



