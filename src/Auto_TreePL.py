#!/usr/bin/env python3
import argparse
import subprocess
import re
import os
import sys
import time

def config_cv(config_file, cv_file, tree_file, num_sites, nthreads, timeline):
    """
    Load the input config file, and generate a config file for cross-validation.
    """
    # Load the input config file and copy the content to the new config file
    with open(config_file, 'r') as src:
        content = src.read()

    with open(cv_file, 'w') as dst:
        dst.write(f"# Your treefile\n")
        dst.write(f"treefile = {tree_file}\n")
        dst.write(f"\n")
        dst.write(f"# smooth value\n")
        dst.write(f"smooth = 100\n")
        dst.write(f"\n")
        dst.write(f"# the number of sites\n")
        dst.write(f"numsites = {num_sites}\n")
        dst.write(f"\n")
        dst.write(f"# random seed\n")
        dst.write(f"seed = 12345\n")
        dst.write(f"\n")
        dst.write(f"# the calibration points\n")
        dst.write(content)
        dst.write(f"\n")
        dst.write(f"# output file name\n")
        dst.write(f"cvoutfile = treepl_cv_results_{timeline}.txt\n")
        dst.write(f"\n")
        dst.write(f"# the number of threads\n")
        dst.write(f"nthreads = {nthreads}\n")
        dst.write(f"\n")
        dst.write(f"# thorough optimization\n")
        dst.write(f"thorough\n")
        dst.write(f"\n")
        dst.write(f"# do cross validation\n")
        dst.write(f"cv\n")
        dst.write(f"# check the calibration\n")
        dst.write(f"checkconstraints\n")

def find_min_chisq_smooth(cv_file):
    """
    Find the optimal smooth value from the cross-validation results.
    """
    min_chisq = None
    corresponding_smooth = None
    
    with open(cv_file, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('chisq:'):
                # retrieve the smooth value and chisq value
                parts = line.split()
                smooth = float(parts[1][1:-1])  # extract the number in the bracket and convert to float
                chisq = float(parts[2])        # extract the chisq value
                
                # check if it is the smallest chisq value
                if min_chisq is None or chisq < min_chisq:
                    min_chisq = chisq
                    corresponding_smooth = smooth
    
    return corresponding_smooth

def config_prime(config_file, prime_file, tree_file, num_sites, nthreads, smooth, timeline):
    """
    Load the input config file, and generate a config file for prime.
    """
    # Load the input config file and copy the content to the new config file
    with open(config_file, 'r') as src:
        content = src.read()

    with open(prime_file, 'w') as dst:
        dst.write(f"# Your treefile\n")
        dst.write(f"treefile = {tree_file}\n")
        dst.write(f"\n")
        dst.write(f"# smooth value\n")
        dst.write(f"smooth = {smooth}\n")
        dst.write(f"\n")
        dst.write(f"# the number of sites\n")
        dst.write(f"numsites = {num_sites}\n")
        dst.write(f"\n")
        dst.write(f"# random seed\n")
        dst.write(f"seed = 12345\n")
        dst.write(f"\n")
        dst.write(f"# the calibration points\n")
        dst.write(content)
        dst.write(f"\n")
        dst.write(f"# output file name\n")
        dst.write(f"primeoutfile = treepl_prime_results_{timeline}.txt\n")
        dst.write(f"\n")
        dst.write(f"# the number of threads\n")
        dst.write(f"nthreads = {nthreads}\n")
        dst.write(f"\n")
        dst.write(f"# thorough optimization\n")
        dst.write(f"thorough\n")
        dst.write(f"\n")
        dst.write(f"# prime mode\n")
        dst.write(f"prime\n")

def run_treepl(config_file, step):
    """Running treePL"""
    try:
        result = subprocess.run(
            ['treePL', config_file],  
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=True  
        )
        return result
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Failed to run TreePL {step} with return code {e.returncode}:")
        print(f"{e.stderr}")
        sys.exit(1)
    except Exception as e:
        print(f"[ERROR] Unexpected error in TreePL {step}: {str(e)}")
        sys.exit(1)

def config_dated(config_file, dated_file, tree_file, num_sites, nthreads, smooth, output_dated_tree, last_five_lines):
    """
    Load the input config file, and generate a config file for cross-validation.
    """
    # Load the input config file and copy the content to the new config file
    with open(config_file, 'r') as src:
        content = src.read()

    with open(dated_file, 'w') as dst:
        dst.write(f"# Your treefile\n")
        dst.write(f"treefile = {tree_file}\n")
        dst.write(f"\n")
        dst.write(f"# smooth value\n")
        dst.write(f"smooth = {smooth}\n")
        dst.write(f"\n")
        dst.write(f"# the number of sites\n")
        dst.write(f"numsites = {num_sites}\n")
        dst.write(f"\n")
        dst.write(f"# random seed\n")
        dst.write(f"seed = 12345\n")
        dst.write(f"\n")
        dst.write(f"# the calibration points\n")
        dst.write(content)
        dst.write(f"\n")
        dst.write(f"# output file name\n")
        dst.write(f"outfile = {output_dated_tree}\n")
        dst.write(f"\n")
        dst.write(f"# the number of threads\n")
        dst.write(f"nthreads = {nthreads}\n")
        dst.write(f"\n")
        dst.write(f"# thorough optimization\n")
        dst.write(f"thorough\n")
        dst.write(f"\n")
        dst.write(f"# do cross validation\n")
        dst.write(f"date\n")
        dst.write(f"log_p\n")
        dst.write(f"\n")
        for line in last_five_lines:
            dst.write(line + "\n")
        dst.write("\n")

def main():
    logo = r"""
     █████╗ ██╗   ██╗████████╗ ██████╗        ████████╗██████╗ ███████╗███████╗██████╗ ██╗
    ██╔══██╗██║   ██║╚══██╔══╝██╔═══██╗       ╚══██╔══╝██╔══██╗██╔════╝██╔════╝██╔══██╗██║
    ███████║██║   ██║   ██║   ██║   ██║ █████╗   ██║   ██████╔╝█████╗  █████╗  ██████╔╝██║
    ██╔══██║██║   ██║   ██║   ██║   ██║ ╚════╝   ██║   ██╔══██╗██╔══╝  ██╔══╝  ██╔═══╝ ██║
    ██║  ██║╚██████╔╝   ██║   ╚██████╔╝          ██║   ██║  ██║███████╗███████╗██║     ██████╗
    ╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝           ╚═╝   ╚═╝  ╚═╝╚══════╝╚══════╝╚═╝     ╚═════╝
                            Script for running TreePL automatically
    """
    print("\033[38;5;46m" + logo + "\033[0m")
    # Initialize argument parser with description
    parser = argparse.ArgumentParser(description="Automated TreePL workflow: Input tree file and calibration config, output time-calibrated tree.")
    
    # Define required command line arguments:
    # -it/--input_treefile: Input tree file in Newick/NEXUS format (required)
    # -ci/--calibration_info: Configuration file with calibration points (required)
    # -ot/--output_tree: Output file name for dated tree (default: "output_dated.tre")
    parser.add_argument("-it", "--input_tree", required=True, help="Input tree file (Newick/NEXUS format)")
    parser.add_argument("-ci", "--calibration_info", required=True, help="Calibration points config file (.config format)")
    parser.add_argument("-ns", "--num_sites", required=True, type=int, help="the number sites from your supermatrix, which is used to estimated the tree above")
    parser.add_argument("-ot", "--output_tree", default="output_dated.tre", help="Path to the output time-tree file (default: output_dated.tre)")
    parser.add_argument("-nt", "--nthreads", default=1, type=int, help="Number of threads for running TreePL (default: 1)")

    # Parse command line arguments
    args = parser.parse_args()

    # Define the timeline
    timeline = time.strftime("%Y%m%d_%H:%M:%S")
    
    # Step 1: Verify input files exist
    # Check if input tree file exists, raise error if not found
    if not os.path.exists(args.input_tree):
        raise FileNotFoundError(f"Input tree file not found: {args.input_tree}")
    # Check if config file exists, raise error if not found
    if not os.path.exists(args.calibration_info):
        raise FileNotFoundError(f"Calibration config file not found: {args.calibration_info}")
    # Check if the number of sites is an integer and greater than 0
    if not isinstance(args.num_sites, int):
        raise ValueError("The number of sites must be an integer")
    if args.num_sites <= 0:
        raise ValueError("The number of sites must be greater than 0")
    # Check if the number of threads is an integer and greater than 0
    if not isinstance(args.nthreads, int):
        raise ValueError("The number of threads must be an integer")
    if args.nthreads <= 0:
        raise ValueError("The number of threads must be greater than 0")
    print("\n" + "="*50)
    print("Starting TreePL Automated Dating Pipeline".center(50))
    print("="*50 + "\n")
    
    # Step 2: Generate temporary cv config file for TreePL
    print("[STEP 1/7] Generating cross-validation config...")
    temp_cv_config = f"treepl_cv_{timeline}.config"
    config_cv(args.calibration_info, temp_cv_config, args.input_tree, args.num_sites, args.nthreads, timeline)
    print(f"✓ Cross-validation config file saved to {temp_cv_config}")

    # Step 3: Run TreePL cross-validation
    print(f"\n[STEP 2/7] Running cross-validation (it may cost some time, using {args.nthreads} threads)...")
    run_treepl(temp_cv_config, "cross-validation")
    print("✓ Cross-validation completed successfully")

    # Step 4: Find the optimal smooth value
    print("\n[STEP 3/7] Finding optimal smooth value...")
    temp_cv_results = f"treepl_cv_results_{timeline}.txt"
    if not os.path.exists(temp_cv_results):
        print(f"✗ Error: Results file {temp_cv_results} not found!")
        sys.exit(1)
    smooth = find_min_chisq_smooth(temp_cv_results)
    print(f"✓ Optimal smooth value found: {smooth}")

    # Step 5: Generate temporary prime config file for TreePL
    print("\n[STEP 4/7] Generating temporary prime config file...")
    temp_prime_config = f"treepl_prime_{timeline}.config"
    config_prime(args.calibration_info, temp_prime_config, args.input_tree, args.num_sites, args.nthreads, smooth, timeline)
    print(f"✓ Prime config file saved to {temp_prime_config}")

    # Step 6: Run TreePL prime
    print(f"\n[STEP 5/7] Running TreePL prime... (it may cost some time, using {args.nthreads} threads)")
    result = run_treepl(temp_prime_config, "prime")
    print("✓ TreePL prime completed successfully")
    output_lines = result.stdout.splitlines()
    last_five_lines = output_lines[-6:] if len(output_lines) >= 6 else output_lines
    
    # Step 7: Generating temporary date config file
    print(f"\n[STEP 6/7] Generating temporary date config file...")
    temp_dated_config = f"treepl_dated_{timeline}.config"
    config_dated(args.calibration_info, temp_dated_config, args.input_tree, args.num_sites, args.nthreads, smooth, args.output_tree, last_five_lines)
    print(f"✓ Dating config file saved to {temp_dated_config}")

    # Step 8: run dating
    print(f"\n[STEP 7/7] Generating the final time tree...")
    run_treepl(temp_dated_config, "dating")
    print(f"✓ Time tree has been saved in {args.output_tree}")

if __name__ == "__main__":
    main()
