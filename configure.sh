#!/bin/bash
singularity build container/bcftools.sif container/bcftools.def 
singularity build container/pypy.sif container/pypy.def 
singularity build container/python3.sif container/python3.def 
singularity build container/shapeit.sif container/shapeit.def 
singularity build container/R.sif container/R.def 
