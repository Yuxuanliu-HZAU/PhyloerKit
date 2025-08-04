#!/usr/bin/env python3
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
import argparse
import subprocess
import re
import os
from Bio import Phylo
from Bio.Phylo import NewickIO
from io import StringIO
from tqdm import tqdm
import tempfile
import shutil
import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from sklearn.linear_model import LinearRegression
from datetime import datetime
import rpy2.robjects as ro
from rpy2.robjects.packages import importr, isinstalled
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter

def print_with_timestamp(message):
    """Print message with timestamp"""
    # If the message is empty or only contains newlines, print an empty line
    if not message.strip():
        print()
    else:
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"[{timestamp}] {message}")

def convert_single_fasta_to_phylip(fasta_file, output_file):
    """
    Convert a single FASTA file to PHYLIP format:
    - Supports species names of any length (Extended PHYLIP format)
    - Ensures ≥2 spaces between species names and sequences
    - Strict length validation
    """
    sequences = {}
    
    # Read FASTA file
    with open(fasta_file, 'r') as f:
        current_id = None
        current_seq = []
        for line in f:
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_id:
            sequences[current_id] = ''.join(current_seq)

    if not sequences:
        raise ValueError(f"No valid sequences found in file: {fasta_file}")

    # Get sequence parameters
    seq_len = len(next(iter(sequences.values())))
    max_name_len = max(len(name) for name in sequences.keys())
    name_width = max(10, max_name_len + 2)  # Minimum 10 characters width

    # Write in PHYLIP format
    with open(output_file, 'w') as out:
        # First line: number_of_sequences sequence_length
        out.write(f"{len(sequences)} {seq_len}\n")
        # Write sequences with proper spacing
        for seq_id, seq in sequences.items():
            out.write(f"{seq_id.ljust(name_width)}{seq}\n")

def process_and_merge_fastas(input_dir, output_file):
    """
    Process all FASTA files in a directory and merge them into a single PHYLIP file
    Steps:
    1. Convert each FASTA to PHYLIP format
    2. Merge all PHYLIP files into one
    """
    # Create temporary directory for intermediate files
    temp_dir = tempfile.mkdtemp()
    try:
        # Step 1: Collect all FASTA files
        fasta_files = [f for f in os.listdir(input_dir) if f.endswith(('.fasta', '.fa'))]
        
        if not fasta_files:
            raise ValueError("No .fasta/.fa files found in the input directory")

        print_with_timestamp(f"Found {len(fasta_files)} FASTA files to process...")
        
        # Step 2: Convert each FASTA to PHYLIP
        phylip_files = []
        # Add timestamp to progress bar description
        timestamp = datetime.now().strftime('[%Y-%m-%d %H:%M:%S]')
        with tqdm(total=len(fasta_files), desc=f"{timestamp} Converting FASTA to PHYLIP") as pbar:
            for fasta_file in fasta_files:
                input_path = os.path.join(input_dir, fasta_file)
                output_path = os.path.join(temp_dir, f"{os.path.splitext(fasta_file)[0]}.phy")
                convert_single_fasta_to_phylip(input_path, output_path)
                phylip_files.append(output_path)
                pbar.update(1)

        # Step 3: Merge all PHYLIP files
        print_with_timestamp("Merging PHYLIP files...")
        with open(output_file, 'w') as outfile:
            for phy_file in phylip_files:
                with open(phy_file, 'r') as infile:
                    # Copy entire content including headers
                    shutil.copyfileobj(infile, outfile)

        print_with_timestamp(f"Conversion complete! Results saved to: {output_file}")
        print_with_timestamp(f"Summary: Processed {len(fasta_files)} FASTA files and merged them into one PHYLIP file")

    finally:
        # Cleanup: Remove temporary files
        shutil.rmtree(temp_dir)

def write_ctl_file(prefix, ndata, seqtype, clock, RootAge, burn_in, sample_freq, num_samples, model, cleandata, print):
    """
    Write the ctl file for running MCMCTree.
    """
    with open("mcmctree_0.ctl", 'w') as ctl_file:
        ctl_file.write(f"          seed = -1\n")
        ctl_file.write(f"       seqfile = {prefix}.phy\n")
        ctl_file.write(f"      treefile = {prefix}.calibrated.tre\n")
        ctl_file.write(f"      mcmcfile = {prefix}.mcmc0.txt\n")
        ctl_file.write(f"       outfile = {prefix}.output0.txt\n\n")
        ctl_file.write(f"         ndata = {ndata}\n")
        ctl_file.write(f"       seqtype = {seqtype}\n")
        ctl_file.write(f"       usedata = 3\n")
        ctl_file.write(f"         clock = {clock}\n")
        ctl_file.write(f"       RootAge = {RootAge}\n\n")
        ctl_file.write(f"         model = {model}\n")
        ctl_file.write(f"         alpha = 0.5\n")
        ctl_file.write(f"         ncatG = 5\n\n")
        ctl_file.write(f"     cleandata = {cleandata}\n\n")
        ctl_file.write(f"       BDparas = 1 1 0.1\n")
        ctl_file.write(f"   kappa_gamma = 6 2\n")
        ctl_file.write(f"   alpha_gamma = 1 1\n\n")
        ctl_file.write(f"   rgene_gamma = 2 20 1\n")
        ctl_file.write(f"  sigma2_gamma = 1 10 1\n\n")
        ctl_file.write(f"      finetune = 1: .1 .1 .1 .1 .1 .1\n\n")
        ctl_file.write(f"         print = {print}\n")
        ctl_file.write(f"        burnin = {burn_in}\n")
        ctl_file.write(f"      sampfreq = {sample_freq}\n")
        ctl_file.write(f"       nsample = {num_samples}\n")

def remove_branch_lengths_and_support(input_tree, output_tree=None):
    """
    Remove branch lengths and support values from a Newick tree while preserving species names in quotes
    
    Parameters:
        input_tree: Input tree file path or Newick string
        output_tree: Output tree file path (optional), if None returns string
    
    Returns:
        If output_tree is not specified, returns the processed Newick string
        If output_tree is specified, returns None (result written to file)
    """
    # Read input tree
    if isinstance(input_tree, str) and '\n' not in input_tree and '(' in input_tree:
        tree_str = input_tree
    else:
        try:
            tree = Phylo.read(input_tree, 'newick')
            tree_str = tree.format('newick')
        except Exception as e:
            raise ValueError(f"Failed to read tree file: {str(e)}")

    # use more precise regular expression pattern matching
    # 1. remove support values (numbers after right parenthesis)
    tree_str = re.sub(r'\)(\d+\.?\d*)', ')', tree_str)
    
    # 2. remove branch lengths (numbers after colon, including scientific notation)
    tree_str = re.sub(r':[0-9.eE+-]+', '', tree_str)
    
    # 3. clean up possible extra spaces, but keep basic structure
    tree_str = re.sub(r'\s+', '', tree_str)
    
    # 4. verify that the tree's parentheses match
    if tree_str.count('(') != tree_str.count(')'):
        raise ValueError("Unmatched parentheses in processed tree")
    
    # 5. ensure the tree ends with a semicolon
    if not tree_str.endswith(';'):
        tree_str += ';'
    
    # 6. verify that the processed tree is still a valid Newick format
    try:
        Phylo.read(StringIO(tree_str), 'newick')
    except Exception as e:
        raise ValueError(f"Invalid Newick format after processing: {str(e)}")

    if output_tree is None:
        return tree_str
    else:
        with open(output_tree, 'w') as f:
            f.write(tree_str)
        return None

def apply_calibrations_to_tree(config_file, tree_file, output_tree, precision, tree_format='newick'):
    calibrations = []
    
    with open(config_file, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        if line.startswith('mrca'):
            # more robust parsing method (handle possible format changes)
            try:
                min_age = None
                max_age = None

                parts = re.split(r'\s*=\s*|\s+', line)  # handle = and space separation
                species1 = parts[-2]
                species2 = parts[-1]
                
                if i + 1 < len(lines):
                    min_line = lines[i+1].strip()
                else:
                    min_line = None

                if i + 2 < len(lines):
                    max_line = lines[i+2].strip()
                else:
                    max_line = None
                
                def extract_age(line):
                    if match := re.search(r'[\d.]+(?:[eE][+-]?\d+)?$', line.strip()):
                        try:
                            return float(match.group())
                        except ValueError:
                            raise ValueError(f"Invalid numerical format: {match.group()} in line: '{line}'")
                    raise ValueError(f"No valid age value found: '{line}'")

                if min_line == None and max_line == None:
                    raise ValueError(f"At least one of the min and max age for the mrca of {species1} and {species2} in the config file should be provided.")

                if min_line != None:
                    if min_line.startswith('min') or min_line.startswith('max'):
                        if min_line.startswith('min'):
                            min_age = extract_age(min_line)
                        elif min_line.startswith('max'):
                            max_age = extract_age(min_line)
                    else:
                        raise ValueError(f"The line next to the mrca line: {line} in the config file is not valid: {min_line}\nIt should start with 'min' or 'max'.")
                
                if min_line != None and max_line != None:
                    if min_line.startswith('min') and max_line.startswith('min'):
                        raise ValueError(f"Two many 'min' in the lines of the config file: {min_line} and {max_line}")

                if min_line != None and max_line != None:
                    if min_line.startswith('max') and max_line.startswith('max'):
                        raise ValueError(f"Two many 'max' in the lines of the config file: {min_line} and {max_line}")
                
                if max_line != None:
                    if max_line.startswith('min') or max_line.startswith('max'):
                        if max_line.startswith('min'):
                            min_age = extract_age(max_line)
                        elif max_line.startswith('max'):
                            max_age = extract_age(max_line)
                            
                if min_age is not None and max_age is not None:
                    if min_age == max_age:
                        raise ValueError(f"The min age and max age for the mrca of {species1} and {species2} in the config file are the same: {min_age} and {max_age}")
                    
                    if min_age > max_age:
                        raise ValueError(f"The min age for the mrca of {species1} and {species2} in the config file is greater than the max age: {min_age} and {max_age}")
                
                if (min_age is not None and min_age < 0) or (max_age is not None and max_age < 0):
                    raise ValueError(f"The min age or max age for the mrca of {species1} and {species2} in the config file is less than 0: {min_age} and {max_age}")

                if min_age is not None:
                    value_minage = round(min_age / 100, precision)

                if max_age is not None:
                    value_maxage = round(max_age / 100, precision)
                
                if min_age is not None and max_age is not None:
                    mrca_value = f" '>{value_minage}<{value_maxage}'".replace('0.', '.')
                elif min_age is not None:
                    mrca_value = f" '>{value_minage}'".replace('0.', '.')
                elif max_age is not None:
                    mrca_value = f" '<{value_maxage}'".replace('0.', '.')
                else:
                    raise ValueError(f"At least the min or max age for the mrca of {species1} and {species2} in the config file should be provided.")
                
                calibrations.append({
                    'species1': species1,
                    'species2': species2,
                    'mrca_value': mrca_value
                })
                
                i += 3
            except (IndexError, ValueError, AttributeError) as e:
                raise ValueError(f"Failed to parse the calibration point, line {i+1}: {str(e)}")
        else:
            i += 1
    
    try:
        tree = Phylo.read(tree_file, tree_format)
    except:
        raise ValueError(f"Failed to read the tree file: {tree_file}")
    
    species_count = len(tree.get_terminals())
    if species_count < 5:
        raise ValueError(f"The tree file {tree_file} is not a valid tree file, the number of species is less than 5")

    def write_annotated(tree, filename, species_count):
        """Process node annotations writing function and add species count in first line"""
        def node_to_str(clade):
            parts = []
            if clade.clades:
                parts.append("(" + ",".join(node_to_str(c) for c in clade.clades) + ")")
            if clade.name:
                parts.append(clade.name)
            if hasattr(clade, 'mrca_value'):
                parts.append(f"{clade.mrca_value}")
            if hasattr(clade, 'branch_length') and clade.branch_length is not None:
                parts.append(f":{clade.branch_length}")
            return "".join(parts)
        
        with open(filename, 'w') as f:
            # Write species count in first line
            f.write(f"{species_count} 1\n")  # Note: using curly braces
            # Write tree structure in second line
            f.write(node_to_str(tree.root) + ";")

    for calibration in calibrations:
        # more efficient terminal node search (avoid traversing the entire tree)
        terminals = []
        for name in [calibration['species1'], calibration['species2']]:
            found = list(tree.find_elements(name=name))
            if not found:
                raise ValueError(f"Failed to find the species '{name}' in the tree file: {tree_file}")
            terminals.extend(found)
        
        if len(terminals) != 2:
            raise ValueError(f"The species pairing is not unique: {calibration['species1']} and {calibration['species2']}")
            
        try:
            mrca = tree.common_ancestor(*terminals)
            setattr(mrca, 'mrca_value', calibration['mrca_value'])  # add attribute
        except Exception as e:
            raise ValueError(f"Failed to find the MRCA node: {calibration}")
            
    write_annotated(tree, output_tree, species_count)
    return tree, species_count

def extract_and_combine_data(file1_path, file2_path, output_path):
    """
    Extract data from two input files and combine them into an output file
    
    Parameters:
        file1_path: the path of the first input file
        file2_path: the path of the second input file
        output_path: the path of the output file
    """
    # extract the first and second columns of the lines starting with "t_n" from the first file
    def extract_from_file1(filepath):
        data = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
            start_processing = False
            for line in lines:
                # check if the line is the start of the data lines
                if "Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width" in line:
                    start_processing = True
                    continue
                
                if start_processing and line.strip() and line.startswith("t_n"):
                    # use regex to split the line, handle multiple spaces
                    parts = re.split(r'\s+', line.strip())
                    if len(parts) >= 2:
                        data.append((parts[0], parts[1]))  # (the first column, the second column)
        return data
    
    # extract the second column of the lines starting with "t_n" from the second file
    def extract_from_file2(filepath):
        data = []
        with open(filepath, 'r') as f:
            lines = f.readlines()
            start_processing = False
            for line in lines:
                if "Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width" in line:
                    start_processing = True
                    continue
                
                if start_processing and line.strip() and line.startswith("t_n"):
                    parts = re.split(r'\s+', line.strip())
                    if len(parts) >= 2:
                        data.append(parts[1])  # the second column
        return data
    
    # extract the data from the two files
    file1_data = extract_from_file1(file1_path)
    file2_data = extract_from_file2(file2_path)
    
    # check if the length of the data is consistent
    if len(file1_data) != len(file2_data):
        raise ValueError(f"The number of t_n entries in the two files {file1_path} and {file2_path} is not consistent")
    
    # write the output file
    with open(output_path, 'w') as out:
        # write the title line
        out.write("X\trun_01\trun_02\n")
        
        # write the data lines
        for (col1, col2), col3 in zip(file1_data, file2_data):
            out.write(f"{col1}\t{col2}\t{col3}\n")

def plot_mcmc_convergence(csv_path, burn_in, sample_freq, num_sample, output_path=None, colormap='coolwarm', tick_step=5):
    """
    Plot MCMC convergence scatter plot
    
    Parameters:
        csv_path: Path to CSV file
        burn_in: Burn-in value
        sample_freq: Sample frequency
        num_sample: Number of samples
        output_path: Path to save the plot (optional)
        colormap: Color scheme for the plot (default: 'coolwarm')
        tick_step: Step size for colorbar ticks (default: 5)
    """
    # 1. Read data
    data = pd.read_csv(csv_path, sep='\t')
    
    # 2. Create figure
    plt.figure(figsize=(8, 8))
    
    # 3. Generate color values from t_n IDs
    node_ids = data['X'].str.extract(r't_n(\d+)').astype(int).values.ravel()
    colors = node_ids
    run1_col = 'run_01'
    run2_col = 'run_02'
    
    # 4. Create scatter plot
    scatter = plt.scatter(data[run1_col], data[run2_col],
                         c=colors,
                         s=120,
                         alpha=0.7,
                         cmap=colormap)
    
    # 5. Set axes
    max_val = max(data[[run1_col, run2_col]].max())
    plt.xlim(0, max_val)
    plt.ylim(0, max_val)
    plt.xlabel('Run 1', fontsize=12)
    plt.ylabel('Run 2', fontsize=12)
    
    # 6. Add two-line title with different font sizes
    ax = plt.gca()
    ax.set_title('MCMCTree Convergence Plot', 
                 fontsize=16, 
                 pad=10,
                 y=1.028,
                 weight='bold')
    
    ax.text(0.5, 1.009,
            f'(burn-in={burn_in}, samplefreq={sample_freq}, nsample={num_sample})',
            horizontalalignment='center',
            verticalalignment='bottom',
            transform=ax.transAxes,
            fontsize=12)
    
    # 7. Add diagonal line
    plt.plot([0, max_val], [0, max_val], 
             color='grey', linewidth=2, linestyle='--')
    
    # 8. Calculate and add R²
    X = data[run1_col].values.reshape(-1, 1)
    y = data[run2_col].values
    model = LinearRegression().fit(X, y)
    r_squared = model.score(X, y)
    # R² position
    text_x = max_val * 0.05
    r2_y = max_val * 0.95
    plt.text(text_x, r2_y, 
             f'R² = {r_squared:.5f}',
             fontsize=12)
    
    # 9. Calculate and add Pearson correlation (below R²)
    corr, p_value = stats.pearsonr(data[run1_col], data[run2_col])
    # r and p-value position
    rp_y = r2_y - max_val * 0.05
    plt.text(text_x, rp_y,
             f'r = {corr:.5f}, p = {p_value:.5e}',
             fontsize=12)
    
    # 10. Add colorbar with t_n ID label and range
    cbar = plt.colorbar(scatter)
    # cbar.set_label('t_n ID', fontsize=12)
    
    # 11. Obtain the minimum and maximum node IDs
    min_id = int(node_ids.min())
    max_id = int(node_ids.max())
    
    # 12. Use the user-specified step size
    tick_positions = np.arange(min_id, max_id + 1, tick_step)
    if tick_positions[-1] != max_id:
        tick_positions = np.append(tick_positions, max_id)
    
    # 13. Create corresponding labels
    tick_labels = [f't_n{id}' for id in tick_positions]
    
    # 14. Set colorbar ticks
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)
    
    # 15. Save or display plot
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    
    # 16. Close figure
    plt.close()
    
    # 17. Print MCMC parameters and convergence status
    print()
    print_with_timestamp("MCMCTree Parameters:")
    print_with_timestamp(f"burn-in: {burn_in}")
    print_with_timestamp(f"samplefreq: {sample_freq}")
    print_with_timestamp(f"nsample: {num_sample}")
    print_with_timestamp("Convergence Analysis:")
    print_with_timestamp(f"R² = {r_squared:.5f}")
    print_with_timestamp(f"r = {corr:.5f}")
    print_with_timestamp(f"p-value = {p_value:.5e}")
    
    if p_value < 0.001 and r_squared >= 0.99:
        print()
        print_with_timestamp("Result: CONVERGED ✓")
        print_with_timestamp("The two MCMC runs show significant convergence.")
    else:
        print()
        print_with_timestamp("Result: NOT CONVERGED ✗")
        print_with_timestamp("Suggestion: Please increase burn-in and nsample values and run this script again to improve convergence.")

def plot_mcmctree(tree_path, output_prefix="mcmctree_plot", 
                 time_correction=100, cex_tips=0.5, lwd_bar=2, 
                 cex_age=0.4, col_age="navy", cex_labels=0.5, 
                 relative_height=0.08, no_margin=True, 
                 label_offset=0.5, formats=('pdf',), ladderize=True): 
    """
    Python wrapper for MCMCtreeR plotting functionality
    
    Parameters:
        tree_path: Path to input tree file
        output_prefix: Output file prefix
        time_correction: Time correction factor
        cex_tips: Size of tip labels
        lwd_bar: Width of uncertainty bars
        cex_age: Size of age labels
        col_age: Color of uncertainty bars
        cex_labels: Size of time axis labels
        relative_height: Relative height of time axis
        no_margin: Whether to remove margins
        label_offset: Offset of labels
        formats: Output formats (currently only supports 'pdf')
        ladderize: Whether to ladderize the tree
    """
    # Use context manager for pandas and numpy conversion
    with localconverter(ro.default_converter + numpy2ri.converter):
        # Check and install required R packages
        required_packages = ['MCMCtreeR', 'ape']
        utils = importr('utils')
        
        for package in required_packages:
            if not isinstalled(package):
                utils.install_packages(package)
        
        # Import required R packages
        mcmctreer = importr('MCMCtreeR')
        grdevices = importr('grDevices')
        base = importr('base')
        
        try:
            # Use R code for reading and plotting
            r_code = f"""
            library(MCMCtreeR)
            
            # Read tree file
            tree <- readMCMCtree("{tree_path}", from.file=TRUE)
            
            # Set up PDF device
            pdf(file="{output_prefix}.pdf", 
                onefile=TRUE, paper="special", useDingbats=FALSE)
            
            # Plot tree
            MCMC.tree.plot(
                phy = tree,
                analysis.type = "MCMCtree",
                time.correction = {time_correction},
                plot.type = "phylogram",
                node.method = "bar",
                cex.tips = {cex_tips},
                lwd.bar = {lwd_bar},
                cex.age = {cex_age},
                col.age = "{col_age}",
                cex.labels = {cex_labels},
                relative.height = {relative_height},
                no.margin = {"TRUE" if no_margin else "FALSE"},
                label.offset = {label_offset},
                add.time.scale = TRUE,
                scale.res = c("Eon", "Period", "Epoch", "Age"),
                ladderize.tree = {"TRUE" if ladderize else "FALSE"}
            )
            
            # Close device
            dev.off()
            """
            
            # Execute R code
            ro.r(r_code)
            print_with_timestamp(f"Generated {output_prefix}.pdf")
                
        except Exception as e:
            print_with_timestamp(f"Error occurred while executing R code: {str(e)}")
            return None

def clean_directory(directory):
    """
    Clean all files in the specified directory
    
    Parameters:
        directory: Path to the directory to be cleaned
    """
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print_with_timestamp(f"Error: Failed to delete {file_path}. Reason: {e}")

def main():
    logo = r"""
     █████╗ ██╗   ██╗████████╗ ██████╗       ███╗   ███╗ ██████╗███╗   ███╗ ██████╗████████╗██████╗ ███████╗███████╗
    ██╔══██╗██║   ██║╚══██╔══╝██╔═══██╗      ████╗ ████║██╔════╝████╗ ████║██╔════╝╚══██╔══╝██╔══██╗██╔════╝██╔════╝
    ███████║██║   ██║   ██║   ██║   ██║█████╗██╔████╔██║██║     ██╔████╔██║██║        ██║   ██████╔╝█████╗  █████╗  
    ██╔══██║██║   ██║   ██║   ██║   ██║╚════╝██║╚██╔╝██║██║     ██║╚██╔╝██║██║        ██║   ██╔══██╗██╔══╝  ██╔══╝  
    ██║  ██║╚██████╔╝   ██║   ╚██████╔╝      ██║ ╚═╝ ██║╚██████╗██║ ╚═╝ ██║ ██████╗   ██║   ██║  ██║███████╗███████╗
    ╚═╝  ╚═╝ ╚═════╝    ╚═╝    ╚═════╝       ╚═╝     ╚═╝ ╚═════╝╚═╝     ╚═╝ ╚═════╝   ╚═╝   ╚═╝  ╚═╝╚══════╝╚══════╝
                            Script for running MCMCTree automatically

    Author: Yuxuan Liu
    Email: 1281096224@qq.com
    Date: 2025-06-04
    Version: 1.0.0
    Description: This script is used to run MCMCTree automatically.
    """
    print()
    print("\033[38;5;39m" + logo + "\033[0m")
    
    parser = argparse.ArgumentParser(description=f"Automated MCMCTree workflow: Takes a phylogenetic tree (Newick format), calibration config file, and a directory of sequence alignments (FASTA format) as input.\n Automatically runs MCMCTree to perform divergence time estimation, generates a time-calibrated tree, and produces a convergence analysis plot to assess MCMC chain reliability.")
    parser.add_argument("-it", "--input_tree", required=True, type=str, help="Input tree file (only Newick format)")
    parser.add_argument("-ic", "--input_config", required=True, type=str, help="Calibration points config file (.config format)")
    parser.add_argument("-o", "--output_dir", default=".", type=str, help="Output directory (default: .)")
    parser.add_argument("-fd", "--fasta_dir", required=True, type=str, help="The directory containing all alignments files (fasta format)")
    parser.add_argument("-p", "--prefix", default="mcmctree", type=str, help="Prefix for the output file (default: mcmctree)")
    parser.add_argument("-bi", "--burn_in", default=200000, type=int, help="Number of burnin for running MCMCTree (default: 200000)")
    parser.add_argument("-sf", "--sample_freq", default=10, type=int, help="Number of sample frequency for running MCMCTree (default: 10)")
    parser.add_argument("-ns", "--num_samples", default=50000, type=int, help="Number of samples for running MCMCTree (default: 50000)")
    parser.add_argument("-rs", "--recommended_sampling", choices=["small", "medium", "high", "super"], help="The recommended sampling scheme for running MCMCTree.")
    parser.add_argument("-ap", "--age_precision", default=4, type=int, help="Precision for the age values in the config file after dividing by 100 (default: 4)")
    parser.add_argument("-ndata", default=5, type=int, help="Number of data for running MCMCTree (default: 5)")
    parser.add_argument("-seqtype", default=0, type=int, help="Sequence type for running MCMCTree (0: nucleotides; 1:codons; 2:AAs; default: 0)")
    parser.add_argument("-clock", default=2, type=int, help="Use clock for running MCMCTree (1: global clock; 2: independent rates; 3: correlated rates default: 2)")
    parser.add_argument("-RootAge", default="'<1'", type=str, help="Age constraint on root age (used if no fossil for root, default: <1)")
    parser.add_argument("-model", default=4, type=int, help="Model for running MCMCTree (0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85; default: 4)")
    parser.add_argument("-cleandata", default=0, type=int, help="Remove sites with ambiguity data (1:yes, 0:no; default: 0)")
    parser.add_argument("-print", default=1, type=int, help="Print the output (0: no mcmc sample; 1: everything except branch rates; 2: everything; default: 1)")
    parser.add_argument("-c", "--colormap", default='coolwarm', type=str, 
                       choices=['viridis', 'plasma', 'inferno', 'magma', 'cividis',
                               'coolwarm', 'rainbow', 'jet', 'turbo'],
                       help="Color scheme for the plot (default: viridis)")
    parser.add_argument("-t", "--tick_step", default=5, type=int,
                       help="Step size for colorbar ticks (default: 5)")

    # Tree display
    parser.add_argument("--plot_tree", action="store_true", default=False,
                       help="Plot the tree (default: False)")
    parser.add_argument("--no_ladderize_tree", action="store_true", default=False,
                       help="Do not ladderize the tree when plotting the MCMCTree timetree (default: False)")
    args = parser.parse_args()

    try:
        # Step 0: change the working directory to the output directory
        print()
        print_with_timestamp("Step 0: Preparing output directory...")
        
        # Create output directory if it doesn't exist
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
            print_with_timestamp("Created output directory")
        else:
            # Clean directory if it exists and is not empty
            if os.listdir(args.output_dir):
                print_with_timestamp("Cleaning output directory...")
                clean_directory(args.output_dir)
                print_with_timestamp("Output directory cleaned")
        
        # Change to output directory
        os.chdir(args.output_dir)
        print_with_timestamp("✓ Working directory prepared successfully")

        # Step 1: Convert the fasta files to the merged phylip file
        print()
        print_with_timestamp("Step 1: Converting FASTA files to PHYLIP format...")
        process_and_merge_fastas(args.fasta_dir, args.prefix + ".phy")
        if not os.path.exists(args.prefix + ".phy"):
            raise RuntimeError("PHYLIP file generation failed")
        print_with_timestamp("✓ FASTA to PHYLIP conversion completed")

        # Step 2: remove the branch lengths and support from the tree
        print()
        print_with_timestamp("Step 2: Processing input tree...")
        tree_str = remove_branch_lengths_and_support(args.input_tree, args.prefix + ".uncalibrated.tre")
        if not os.path.exists(args.prefix + ".uncalibrated.tre"):
            raise RuntimeError("Tree processing failed")
        print_with_timestamp("✓ Tree processing completed")

        # Step 3: Calibrate the tree automatically
        print()
        print_with_timestamp("Step 3: Performing automatic tree calibration...")
        try:
            with open(args.input_config, 'r') as config_file:
                config_content = config_file.read()
            config_content = config_content.replace("'", "")
            new_config_file_path = args.prefix + ".calibrated.config"
            with open(new_config_file_path, 'w') as new_config_file:
                new_config_file.write(config_content)
            args.input_config = new_config_file_path
            
            tree_str = remove_branch_lengths_and_support(args.input_tree, args.prefix + ".uncalibrated.tre")
            tree, species_count = apply_calibrations_to_tree(
                config_file=args.input_config,
                tree_file=args.prefix + ".uncalibrated.tre",
                output_tree=args.prefix + ".calibrated.tre",
                precision=int(args.age_precision)
            )
            print_with_timestamp("✓ Automatic tree calibration completed")
            print_with_timestamp(f"The calibrated tree is saved in the file: {os.path.join(args.output_dir, f'{args.prefix}.calibrated.tre')}")
        except Exception as e:
            raise RuntimeError(f"Automatic tree calibration failed: {str(e)}")

        # Step 4: write the ctl file for running MCMCTree
        print()
        print_with_timestamp("Step 4: Writing the control file for MCMCTree (0th round)...")
        write_ctl_file(args.prefix, args.ndata, args.seqtype, args.clock, args.RootAge,
                      args.burn_in, args.sample_freq, args.num_samples, args.model,
                      args.cleandata, args.print)
        print_with_timestamp("✓ Control file for MCMCTree (0th round) written successfully")

        # Step 5: Run MCMCTree (preparation) to generate the out.BV file
        print()
        print_with_timestamp("Step 5: Running MCMCTree (preparation) to generate the out.BV file")
        result = subprocess.run(["mcmctree", "mcmctree_0.ctl"], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"MCMCTree preparation running failed: {result.stderr}")
        print_with_timestamp("✓ MCMCTree preparation running completed")

        # Step 6: Prepare and run MCMCTree (1st round)
        print()
        print_with_timestamp("Step 6: Running MCMCTree (1st round)...")
        try:
            if os.path.exists("out.BV"):
                os.rename("out.BV", "in.BV")
                shutil.copy("in.BV", "out.BV")
                shutil.copy("mcmctree_0.ctl", "mcmctree_1.ctl")

                with open("mcmctree_1.ctl", 'r') as f:
                    content = f.read()
                content = content.replace("usedata = 3", "usedata = 2 * use in.BV directly")
                content = content.replace(f"outfile = {args.prefix}.output0.txt", 
                                       f"outfile = {args.prefix}.output1.txt * use in.BV directly")
                content = content.replace(f"mcmcfile = {args.prefix}.mcmc0.txt", 
                                       f"mcmcfile = {args.prefix}.mcmc1.txt * use in.BV directly")
                with open("mcmctree_1.ctl", 'w') as f:
                    f.write(content)
                
                result = subprocess.run(["mcmctree", "mcmctree_1.ctl"], capture_output=True, text=True)
                if result.returncode != 0:
                    raise RuntimeError(f"MCMCTree 2nd round failed: {result.stderr}")
                print_with_timestamp("✓ MCMCTree 1st round completed")
            else:
                raise RuntimeError("out.BV file not found after first MCMCTree run")
        except Exception as e:
            raise RuntimeError(f"MCMCTree 1st round preparation or execution failed: {str(e)}")

        # Step 7: run MCMCTree (2nd round)
        print()
        print_with_timestamp("Step 7: Running MCMCTree (2nd round)...")
        try:    
            shutil.copy("mcmctree_1.ctl", "mcmctree_2.ctl")
            with open("mcmctree_2.ctl", 'r') as f:
                content = f.read()
            content = content.replace(f"outfile = {args.prefix}.output1.txt", 
                                   f"outfile = {args.prefix}.output2.txt * use in.BV directly")
            content = content.replace(f"mcmcfile = {args.prefix}.mcmc1.txt", 
                                   f"mcmcfile = {args.prefix}.mcmc2.txt * use in.BV directly")
            with open("mcmctree_2.ctl", 'w') as f:
                f.write(content)
            result = subprocess.run(["mcmctree", "mcmctree_2.ctl"], capture_output=True, text=True)
            if result.returncode != 0:
                raise RuntimeError(f"MCMCTree 3rd round failed: {result.stderr}")
            if not os.path.exists("FigTree.tre"):
                raise RuntimeError("The time tree was not normally generated during MCMCTree execution")
            os.rename("FigTree.tre", f"{args.prefix}.MCMCTree_dated.tre")
            print_with_timestamp("✓ MCMCTree 3rd round completed")
        except Exception as e:
            raise RuntimeError(f"MCMCTree 3rd round preparation or execution failed: {str(e)}")

        # Step 8: Extract the data from the output files for convergence analysis
        print()
        print_with_timestamp("Step 8: Extracting and analyzing output data...")
        if not (os.path.exists(args.prefix + ".output1.txt") and 
                os.path.exists(args.prefix + ".output2.txt")):
            raise RuntimeError("Required output files not found")
        extract_and_combine_data(args.prefix + ".output1.txt", 
                               args.prefix + ".output2.txt",
                               args.prefix + ".convergence_analysis.tsv")
        print_with_timestamp("✓ Data extraction completed")

        # Step 9: Execute the convergence analysis and plot the results
        print()
        print_with_timestamp("Step 9: Executing the convergence analysis and plot the results...")       
        plot_mcmc_convergence(args.prefix + ".convergence_analysis.tsv",
                              args.burn_in, args.sample_freq, args.num_samples,
                              args.prefix + ".convergence_analysis.png",
                              args.colormap, args.tick_step)
        if not os.path.exists(args.prefix + ".convergence_analysis.png"):
            raise RuntimeError(f"Failed to execute the convergence analysis and plot the results.")
        else:
            print_with_timestamp(f"✓ Convergence analysis execution and plotting completed.\n  The results are plotted in {os.path.join(args.output_dir, f'{args.prefix}.convergence_analysis.png')}")
            tmp_files = [f for f in os.listdir('.') if f.startswith('tmp')]
            for f in tmp_files:
                if os.path.exists(f):
                    os.remove(f)

        print()
        print_with_timestamp("✓ All steps completed successfully!")

        # Step 10: Plot the tree
        print()
        print_with_timestamp("Step 10: Plotting the tree...")
        if args.plot_tree:
            try:
                tree = Phylo.read(args.prefix + ".MCMCTree_dated.tre", "nexus")
                tree_path = args.prefix + ".MCMCTree_dated.tre"
                if args.no_ladderize_tree:
                    plot_mcmctree(tree_path=tree_path,
                                output_prefix=args.prefix + ".MCMCTree_dated",
                                ladderize=False)
                else:
                    plot_mcmctree(tree_path=tree_path,
                                output_prefix=args.prefix + ".MCMCTree_dated",
                                ladderize=True)
                
                if not os.path.exists(args.prefix + ".MCMCTree_dated.pdf"):
                    raise RuntimeError(f"\n❌ Tree plotting failed")
                else:
                    print_with_timestamp(f"✓ Tree plotting completed.\n  The tree is plotted in {os.path.join(args.output_dir, f'{args.prefix}.MCMCTree_dated.pdf')}")
            except Exception as e:
                raise RuntimeError(f"Tree plotting failed: {str(e)}")
        else:
            print_with_timestamp(f"The tree is not plotted.\n  Please run the script with the --plot_tree option to plot the tree.")

    except Exception as e:
        print()
        print_with_timestamp(f"❌ Error: {str(e)}")
        print_with_timestamp("Pipeline execution failed. Please check the error message above.")
        tmp_files = [f for f in os.listdir('.') if f.startswith('tmp')]
        for f in tmp_files:
            if os.path.exists(f):
                os.remove(f)
        sys.exit(1)

if __name__ == "__main__":
    main()
