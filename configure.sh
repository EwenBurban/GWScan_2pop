#!/bin/bash
sudo singularity build container/bcftools.sif container/bcftools.def 
sudo singularity build container/pypy.sif container/pypy.def 
sudo singularity build container/python3.sif container/python3.def 
sudo singularity build container/shapeit.sif container/shapeit.def 
sudo singularity build container/R.sif container/R.def 
