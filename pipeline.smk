
## import module and info from config file
import os
import sys
import re
wdir=config['work_dir']
binpath=config['binpath']
Nthread_shape_it=8
popfile=config['popfile']
locusLength=int(config['locusLength'])
nameA=config['nameA']
nameB=config['nameB']
contig_file=config['contig_file']
## list vcf files ##
vcf_list_files = os.listdir(wdir)
pattern= re.compile('.vcf')
vcf_list_files = [ x for x in vcf_list_files if pattern.search(x) ]
## Singularity config ##
container_path =binpath + '/container' 
Sc='singularity exec --bind {0},{1} {2}'.format(binpath,wdir,container_path)

## pipeline ##
wildcard_constraints:
    wdir=wdir,
    vcf='|'.join(vcf_list_files)

rule targets:
    input:
        rawdata = expand('{wdir}/{vcf}',wdir=wdir,vcf=vcf_list_files),
        abc_stats = expand('{wdir}/ABCstat_locus.txt',wdir=wdir)
    shell:
        """
        rm shapeit*
        """

rule format_vcf:
    input:
        '{wdir}/{vcf}'
    output:
        temp('{wdir}/formated/{vcf}')
    shell:
        """
        {Sc}/bcftools.sif bcftools annotate -x 'FORMAT' {input} | cat > {output}        
        """

rule filtering:
    input:
        '{wdir}/formated/{vcf}'
    output:
        temp('{wdir}/filtered/{vcf}')
    shell:
        """
        {Sc}/bcftools.sif bcftools view --max-alleles 2 --exclude-types indels {input} | cat > {output}
        """

rule phasing:
    input:
        '{wdir}/filtered/{vcf}'
    output:
        temp('{wdir}/phased/{vcf}')
    threads: Nthread_shape_it
    shell:
        """
        {Sc}/shapeit.sif shapeit --input-vcf {input} -O {output} -T {Nthread_shape_it} 
        {Sc}/shapeit.sif shapeit -convert --input-haps {output} --output-vcf {output}
        """

rule aggregate:
    input: expand('{wdir}/phased/{vcf}',wdir=wdir,vcf=vcf_list_files)
    output: '{wdir}/all_phased.vcf'
    shell:
        """
        a=({wdir}/phased/*)
        grep -h '#' $a > {output}
        grep -hv '#' {wdir}/phased/*.vcf >> {output}
        """


rule haplotyping:
    input:
        '{wdir}/all_phased.vcf'
    output:
        '{wdir}/haplotyped.vcf'
    resources:
        mem_mb=2*len(vcf_list_files)*2000
    shell:
        """
        {Sc}/python.sif python3 {binpath}/make_haplo_vcf.py {input} {output}
        """


rule calculate_globalstat:
    input:
        vcf = '{wdir}/haplotyped.vcf',
        popfile = expand('{popfile}',popfile=popfile),
        contig_file = expand('{contig_file}',contig_file=contig_file)
    output:
        '{wdir}/ABCstat_locus.txt'
    shell:
        """
		 {Sc}/python.sif python3 {binpath}/vcf2abc.py data={input.vcf} \
			popfile={input.popfile} contig_file={input.contig_file} nameA={nameA}\
			nameB={nameB} window_size={locusLength} output_dir={wdir}
        """

rule visualisation:
    input:
        '{wdir}/ABCstat_loci.txt'
    output:
        '{wdir}/gw_plot.pdf'
    shell:
        """
        {Sc}/R.sif Rscript {binpath}/do_plot.R input_file={input} output={output}
        """
