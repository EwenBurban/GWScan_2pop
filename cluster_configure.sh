#!/bin/bash
sudo singularity build --remote container/bcftools.sif container/bcftools.def 
sudo singularity build --remote container/python.sif container/python.def 
sudo singularity build --remote container/shapeit.sif container/shapeit.def 
sudo singularity build --remote container/R.sif container/R.def 
