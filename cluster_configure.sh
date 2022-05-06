#!/bin/bash
singularity build --remote container/bcftools.sif container/bcftools.def 
singularity build --remote container/python.sif container/python.def 
singularity build --remote container/shapeit.sif container/shapeit.def 
singularity build --remote container/R.sif container/R.def 
