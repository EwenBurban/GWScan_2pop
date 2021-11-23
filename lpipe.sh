#!/usr/bin/bash
#SBATCH --mail-user=ewen.burban@univ-rennes1.fr
#SBATCH --mail-type=all
#SBATCH --time=12-60:00:00
## launch the pipeline of GWScan_2pop 
## the provided argument is for --configfile, expecting the yaml file

## load the necessary environement (used for Genouest cluster)
. /local/env/envsnakemake-6.0.5.sh
. /local/env/envsingularity-3.8.0.sh

## directory of the code
binpath='/home/genouest/cnrs_umr6553/eburban/GW_DILS/input_module'

## cluster version with slurm | If you use the cluster version, you must edit the filet cluster.json
snakemake --snakefile ${binpath}/pipeline.smk -p -j 35 --configfile ${1} --cluster-config ${binpath}/cluster.json --cluster "sbatch --nodes={cluster.node} --ntasks={cluster.n} --cpus-per-task={cluster.cpusPerTask} --time={cluster.time} --mem-per-cpu={cluster.memPerCpu} -p {cluster.p} "  --latency-wait 60 --nolock 
## without slurm
#snakemake --snakefile ${binpath}/pipeline.smk -p -j 35 --configfile ${1}  --latency-wait 60 --nolock 
