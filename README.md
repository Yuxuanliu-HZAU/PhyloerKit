# PhyloerKit
- [Auto-MCMCTree.py](https://github.com/SunLab-MiaoPu/PhyloerKit#auto-mcmctreepy)

---

## [Auto-MCMCTree.py](https://github.com/SunLab-MiaoPu/PhyloerKit/blob/main/src/Auto_MCMCTree.py)

![logo-auto-mcmc](https://github.com/SunLab-MiaoPu/PhyloerKit/blob/main/images/auto-mcmctree_logo.png)

### [Manual MCMCTree steps](https://github.com/sabifo4/Tutorial_MCMCtree/tree/main) vs `Auto-MCMCTree.py`

| Manual MCMCTree Steps | This Script |
|-----------------------|-------------|
| 1. Manually convert FASTA to PHYLIP manually | ✅ Automatic conversion |
| 2. Manually calibrate the treefile manually(troublesome and time-consuming) | ✅ Automatically calibrate the tree  according to the user-provided calibration config file |
| 3. Manually create 3 control files | ✅ Automatically generate all CTL files for running MCMCTree according to the user-specified parameters|
| 4. Manually run `mcmctree` 1 time for preparation and 2 times for tree inference | ✅ Automatically handle all executions |
| 5. Manually execute convergence analysis | ✅ Automatically calculate R² to judgle the convergence |
| 6. Manually generate figures of convergence analysis and timetree with separate and additional tools | ✅ Automatic visualization |

**Typical manual workflow requires 10+ commands** → **This script reduces to ONE command to execute the full pipeline!!!**
**Thus, try to use this script to run MCMCTree automatically!!!**

## Dependencies

- Python 3.6+
- MCMCTree (PAML package)
- BioPython
- pandas
- matplotlib
- rpy2 (for tree plotting)
- R package "MCMCtreeR"

Install dependencies:
```bash
conda create -n mcmctree python=3.10
conda activate mcmctree
conda install -c conda-forge numpy=1.24.3 -y
conda install -c conda-forge scipy=1.10.1 -y
conda install -c conda-forge pandas=2.0.3 -y
conda install -c conda-forge biopython=1.81 -y
conda install -c conda-forge matplotlib=3.7.1 -y
conda install -c conda-forge tqdm=4.65.0 -y
conda install -c conda-forge scikit-learn=1.3.0 -y
conda install bioconda::paml
pip install rpy2
conda install conda-forge::r-base
R
> install.packages("MCMCtreeR")
> install.packages("ape")
> quit()
```

## Pipeline overview

This Python script automates the complete MCMCTree workflow for Bayesian divergence time estimation. It handles:

- File format conversion (FASTA → PHYLIP)
- Tree calibration with fossil constraints
- MCMCTree execution (1st round for preparation and 2 rounds for timetree inference)
- Convergence diagnostics and results visualization
- Visualization of time-calibrated trees


## Input files

> [!TIP] 
> Only three input files are needed for running this script!

1. **Input tree** (`-it`):
   - Newick format with branch lengths/support values or not
   - Tip names must match FASTA headers and calibration file

2. **Calibration config file** (`-ic`):
   - The format of the calibration config file is showed as follows:
   ```text
   # Calibration point 1
   mrca = <calibration_name> <species1> <species2>
   min = <calibration_name> 6.5 # the minmum calibration age
   max = <calibration_name> 10.0 # the maximum calibration age

   # Calibration point 2
   mrca = <calibration_name> <species1> <species2>
   min = <calibration_name> <age> # the minmum calibration age
   
   # Calibration point 3
   mrca = <calibration_name> <species1> <species2>
   max = <calibration_name> <age> # the minmum calibration age
   ...
   ```

3. **FASTA directory** (`-fd`):
   - Contains one FASTA per locus
   - All sequences must be aligned
   - Example:
     ```
     >Eleagnus_mollifolia
     ATGC...
     >Eleagnus_angustata
     ATGC...
     ```


## Basic usage:
```bash
python Auto-MCMCTree.py \
    -it input.tre \
    -ic calibrations.config \
    -fd fasta_directory/ \
    -o results/
```

## Test run (with the [`example_data`](https://github.com/SunLab-MiaoPu/PhyloerKit/tree/main/example_dataset)):

#### Installation (with the directory `example_dataset`)
```
git clone https://github.com/SunLab-MiaoPu/PhyloerKit
```

#### Step1: Change the working directory
```
cd ./PhyloerKit/
```

#### Step2: Run the script
We prepared two suites of input files (with different config files for calibration)

- Run the first type of input files
```bash
python ./src/Auto_MCMCTree.py \
-it ./example_dataset/Test-Auto-MCMCTree/input_ASTRAL_HybSuite_RLWP_sorted_rr.tre \
-ic ./example_dataset/Test-Auto-MCMCTree/input1_Elaeagnaceae.config \
-fd ./example_dataset/Test-Auto-MCMCTree/input1_alignments/ \
-o ./example_dataset/Test-Auto-MCMCTree/output1/ \
--plot_tree -p Eleagnus_1
```

- Run the second type of input files
```bash
python ./src/Auto_MCMCTree.py \
-it ./example_dataset/Test-Auto-MCMCTree/input_ASTRAL_HybSuite_RLWP_sorted_rr.tre \
-ic ./example_dataset/Test-Auto-MCMCTree/input2_Elaeagnaceae.config \
-fd ./example_dataset/Test-Auto-MCMCTree/input2_alignments/ \
-o ./example_dataset/Test-Auto-MCMCTree/output2/ \
--plot_tree -p Eleagnus_2
```

## 🛠️ Core parameters

### Execution control:
```bash
-bi 200000      # Burn-in (default: 200000)
-sf 10          # Sampling frequency (default: 10)
-ns 50000       # Number of samples (default: 50000)
```

### Model settings:
```bash
-model 4        # Substitution model (0-4, default: 4=HKY85)
-clock 2        # Clock model (1-3, default: 2=independent rates)
-seqtype 0      # 0=nucleotides, 1=codons, 2=AAs
```

## 📊 Output files

The output files of type1 are showed below:
```
results/
├── Elaeagnus_1.phy                    # Merged PHYLIP alignment as input file for running MCMCTree
├── Elaeagnus_1.calibrated.tre         # Calibrated Newick tree as input file for running MCMCTree
├── Elaeagnus_1.MCMCTree_dated.tre     # Final time tree generated via MCMCTree (NEXUS)
├── Elaeagnus_1.convergence_analysis.png  # The figure displaying the covergence analysis results
└── Elaeagnus_1.MCMCTree_dated.pdf     # Time tree visualization
```

### `Elaeagnus_1.convergence_analysis.png`
This figure displays the convergence results:
<div align="center">
  <img src="https://github.com/SunLab-MiaoPu/PhyloerKit/raw/main/images/Eleagnus_1.convergence_analysis.png" 
       alt="Convergence Analysis Plot" 
       style="width:60%; max-width:800px;">
</div>


### `Elaeagnus_1.MCMCTree_dated.pdf`
**This pdf file displays the dated tree file:**
<br>

<div align="center">
  <img src="https://github.com/SunLab-MiaoPu/PhyloerKit/blob/main/images/Eleagnus_1.MCMCTree_dated.png" 
       alt="Dated tree Plot" 
       style="width:60%; max-width:800px;">
</div>


> [!TIP]
> - For precise editing, you may import this file into Adobe Illustrator, adjust the font size or the font type, etc:
<br>

<div align="center">
  <img src="https://github.com/SunLab-MiaoPu/PhyloerKit/blob/main/images/Eleagnus_1.MCMCTree_dated_2.png" 
       alt="Dated tree Plot 2" 
       style="width:60%; max-width:800px;">
</div>

## 💡 Advanced options

### Recommended sampling Schemes:
```bash
-rs small    # Quick test (bi=10000, ns=10000)
-rs medium   # Medium (bi=100000, ns=20000)
-rs high     # Production (bi=200000, ns=50000)
```

### Visualization control:
```bash
--plot_tree          # Generate PDF tree visualization
--no_ladderize_tree  # Disable automatic ladderizing
-c viridis           # Color scheme for plots
-t 10                # Node label frequency on colorbar
```

## ⚠️ Troubleshooting

**Common Issues:**
1. *"Species not found in tree"*:
   - Verify tip names match exactly between tree, FASTA, and config

2. *Low R² in convergence plot*:
   - Increase burn-in and sample size:
     ```bash
     -bi 500000 -ns 100000
     ```

## Full pipeline parameters

```
options:
  -h, --help            show this help message and exit
  -it, --input_tree INPUT_TREE
                        Input tree file (only Newick format)
  -ic, --input_config INPUT_CONFIG
                        Calibration points config file (.config format)
  -o, --output_dir OUTPUT_DIR
                        Output directory (default: .)
  -fd, --fasta_dir FASTA_DIR
                        The directory containing all alignments files (fasta format)
  -p, --prefix PREFIX   Prefix for the output file (default: mcmctree)
  -bi, --burn_in BURN_IN
                        Number of burnin for running MCMCTree (default: 200000)
  -sf, --sample_freq SAMPLE_FREQ
                        Number of sample frequency for running MCMCTree (default: 10)
  -ns, --num_samples NUM_SAMPLES
                        Number of samples for running MCMCTree (default: 50000)
  -rs, --recommended_sampling {small,medium,high,super}
                        The recommended sampling scheme for running MCMCTree.
  -ap, --age_precision AGE_PRECISION
                        Precision for the age values in the config file after dividing by 100 (default: 4)
  -ndata NDATA          Number of data for running MCMCTree (default: 5)
  -seqtype SEQTYPE      Sequence type for running MCMCTree (0: nucleotides; 1:codons; 2:AAs; default: 0)
  -clock CLOCK          Use clock for running MCMCTree (1: global clock; 2: independent rates; 3: correlated rates default: 2)
  -RootAge ROOTAGE      Age constraint on root age (used if no fossil for root, default: <1)
  -model MODEL          Model for running MCMCTree (0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85; default: 4)
  -cleandata CLEANDATA  Remove sites with ambiguity data (1:yes, 0:no; default: 0)
  -print PRINT          Print the output (0: no mcmc sample; 1: everything except branch rates; 2: everything; default: 1)
  -c, --colormap {viridis,plasma,inferno,magma,cividis,coolwarm,rainbow,jet,turbo}
                        Color scheme for the plot (default: viridis)
  -t, --tick_step TICK_STEP
                        Step size for colorbar ticks (default: 5)
  --plot_tree           Plot the tree (default: False)
  --no_ladderize_tree   Do not ladderize the tree when plotting the MCMCTree timetree (default: False)
```
