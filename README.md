# Description
This is a repository for the snakemake version of the bash GWAS pipeline compatible with the Clemson University's Center for Human Genetics (CUCHG) High Performance Computing (HPC) cluster. This version uses kinship matrix calculated from GEMMA. To use GCTA's ```--make-grm-inbred``` flag for GRM, swap ```Snakefile``` for ```Snakefile_GCTA``` in the ```snakemake_submitter.sh```.

![Pipeline Schematic](https://github.com/vshanka23/dgrp_gwas_final/blob/main/pipeline_schematic.jpg)

- ```slurm/config.yaml```: config file for HPC architecture and slurm compatibility
- ```snakemake_submitter.sh```: initiates conda environment and submits the snakemake job to snakemake
- ```initiator.sh```: sets up the directory and launches the *snakemake_submitted.sh*
- ```Snakefile```: the pipeline
- ```config.yaml```: environmental variables for the pipeline
- ```inputs/```: directory with covariates and phenotype data
- ```scripts/```: directory with R scripts for processing, FDR adjustment and figure generation

# Prerequisites 
***only install these if not running the pipeline on CUCHG's HPC***

- [Anaconda3/miniconda3](https://docs.anaconda.com/anaconda/install/linux/)
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [slurm](https://slurm.schedmd.com/sbatch.html)
- [BCFtools](https://samtools.github.io/bcftools/)
- [R](https://www.r-project.org)
- [PLINK 1.9](https://www.cog-genomics.org/plink/)
- [GCTA](https://yanglab.westlake.edu.cn/software/gcta/)
- [GEMMA](https://github.com/genetics-statistics/GEMMA)
- [OpenBLAS](https://www.openblas.net)
- [GSL](https://www.gnu.org/software/gsl/)

# Instructions
## Adding user- and project- specific information
***generally, add information encompassed by "<>" in the files below***
- ```slurm/config.yaml```: 
    - add max number of jobs ([integer](https://en.wikipedia.org/wiki/Integer_(computer_science)))
    - partition name ([string](https://wlm.userweb.mwn.de/Stata/wstavart.htm))
    - max number of cpus (integer) and max RAM (integer in MB): contact systems administrators for these values and do not edit once established
- ```config.yaml```: Add path to working directory for DEST after removing ```<>```.
- ```snakemake_submitter.sh```:
    - sbatch parameters: 
        - add job name (string)
        - partition name (string)
        - modify time if needed (in Hr:Min:Sec format)
        - output and error (add path to working directory, same as ```DEST``` from ```config.yaml```, but leave the ```/log```... parts unchanged)
        - mail-user (add user email address)
    - ```cd``` line: add path to working directory (same as ```DEST``` from ```config.yaml```)
    - ***Only change this if you are not running on Secretariat.*** ```source``` line: add path to conda initiation script (```conda.sh```) to choose the right conda. There is also a mamba init line if you use mamba. Remove if there is no mamba installation.
    - conda activate line: add the name of the environment with a working snakemake installation (on Secretariat it is "snakemake").

## Prepare the phenotype file

The phenotype file(s) must be placed within ```inputs/pheno/```. The phenotype file must be tab-delimited in this format:

```
DGRP-3-0001	DGRP-3-0001	19.5
DGRP-3-0002	DGRP-3-0002	22.9
DGRP-3-0003	DGRP-3-0003	19.9
DGRP-3-0004	DGRP-3-0004	21.1
DGRP-3-0005	DGRP-3-0005	15.8
DGRP-3-0006	DGRP-3-0006	18.2
DGRP-3-0007	DGRP-3-0007	20.4
DGRP-3-0008	DGRP-3-0008	15.8
DGRP-3-0009	DGRP-3-0009	24.5
DGRP-3-0010	DGRP-3-0010	23.6
```
Phenotype files are automatically captured and fed into the workflow.

1. Open ssh shell (using [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html) or [Putty](https://www.putty.org/) on Windows, terminal app on Mac OS and Linux) and connect to the head/master/login node
2. Make a working directory for the analysis and ```git clone``` this repository [(```git clone https://github.com/vshanka23/dgrp_gwas_final.git```)]
3. Copy ```Snakefile```, ```snakemake_submitter.sh```, ```config.yaml```, ```slurm/config.yaml``` and ```initiator.sh``` to working directory
4. Make sure the variables encompassed by ```<>``` in ```slurm/config.yaml```, ```config.yaml``` and ```snakemake_submitter.sh``` have been modified to reflect info specific to your run (eg: working directory, raw data location, etc)
5. Open a ssh shell and run:
    ```
    ##Initialize the correct conda and bring conda into bash environment
    source <path to conda initialization script>
    ##Activate the correct conda environment containing snakemake installation
    conda activate <snakemake conda environment>

    cd <working directory containing analysis pipeline files>
    ```
    Generate DAG figure:
    ```
    snakemake -n -p -s Snakefile --configfile config.yaml --profile slurm --dag | display | dot
    ```
    Generate the workflow:
    ```
    snakemake -n -p -s Snakefile --configfile config.yaml --profile slurm
    ```

## II. Actual run (head/master/login node)

If step 5 in test run (*Generate the DAG figure* and *Generate the workflow* commands) do not generate any errors (red text), run:
```
./initiator.sh
```

# Tracking progress
There are three places to check for progress:
1. ```squeue```
2. This pipeline (when run successfully) will create ```log``` and ```logs_slurm``` directories within the working directory. In the ```log``` directory, look for ```output_<job_ID>.txt``` and ```error_<job_ID>.txt``` for current status of the run. When the run is successful, the last line should contain a ```x of x steps (100%) done```.
3. In the ```logs_slurm``` directory, the most current log files with specific rule names on the file names represent the current statuses.

# Future additions
1. Need to add annotations from VEP on DGRP3.