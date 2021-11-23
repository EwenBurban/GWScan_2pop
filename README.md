# GWScan_2pop
This pipeline aim to do genome wide scan of various measures of  molecular population genomics. 
It suppose that you have vcf files containing 2 populations.
This pipeline will do the following thing : 
1. Formating vcf file in order to keep only allelic information
2. Phasing haplotypes
3. For each individual, randomly sample one haplotype
4. calculate for each window the statistics => ABCstat_loci.txt
5. create a simple manhatan plot for each statistics in => gw_plot.pdf
 
## Prerequisite
The following list of programs must be install on your machine:
- Snakemake
- Singularity

## Installation
```
git clone https://github.com/EwenBurban/GWScan_2pop.git
cd GW_Scan
chmod a+x configure.sh
./configure.sh
```
## Usage
Before launching the pipeline you must enter various information in 2 files :
- the config file (by default I call it config.yaml) 
- the popfile

an exemple of those two file is aviable (exemple_config.yaml and exemple_popfile.csv)

### config file

The config file provide various information that are mandatory for the pipeline

- work_dir: #the complete path dir where your vcf files are stored
- binpath: #the complete path where GW_Scan is installed
- popfile: #the complete path of your the popfile
- locusLength: #window size of analyse
- nameA: #name of one species (must be the same name in popfile)
- nameB: # same here
- Nref: #this argument is a remainder of part of DILS script. You don’t need to put an exact value (by default I put 85000)
- mu: #this argument is a remainder of part of DILS script. You don’t need to put an exact value (by default I put 3e-8)
- rho_over_theta: #this argument is a remainder of part of DILS script. You don’t need to put an exact value (by default I put 1.0)

### popfile
The popfile inform the pipeline about the population attribution of each individual that you choiced to put in your analyse.
It’s wrote using the csv convention.

### Edit lpipe.sh
lpipe.sh is the script that launch the pipeline. Before launching you must :
- edit the binpath value (put the same dir than in the config file)
- put in comment the line that is you wont use depending of  your cluster management system (with or without slurm)
    - NB : if you use slurm in your cluster, you must edit the cluster.json file
### launch the pipeline

without slurm
```
lpipe.sh <config.yaml path>
```
with slurm
```
sbatch lpipe.sh <config.yaml path>
```

