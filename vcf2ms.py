import pandas as pd
import sys
import re

argv={x.split('=')[0]: x.split('=')[1] for x in sys.argv[1:]}

vcf_file=argv['vcf_file']
popfile=argv['popfile']
bpfile_path=argv['bpfile_path']
outputfile_path=argv['outputfile_path']
locusLength=int(argv['locusLength'])
nameA=argv['nameA']
nameB=argv['nameB']
Nref=int(argv['Nref'])
mu=float(argv['mu'])
rho_over_theta=float(argv['rho_over_theta'])

vcf = pd.read_csv(vcf_file,sep='\t',comment='#',header=None)

# grep header and apply it to vcf data
pattern = re.compile('#CHROM')
with open(vcf_file,'r') as f:
    line = f.readline()
    while line : 
        if pattern.search(line) :
            header=line.strip().split('\t')
            line=False
        else:
            line=f.readline()

vcf.set_axis(header,axis=1,inplace=True) # rename columns name to header

# subset wild and domesticated population
pop = pd.read_csv(popfile,sep=',',header=0)
popA = [x for x in set( pop[nameA].to_list()) if str(x) != 'nan']
popB = [x for x in set( pop[nameB].to_list()) if str(x) != 'nan']

wild_dom = popA + popB

wild_dom= wild_dom +['#CHROM','POS'] # add the first two columns that contains 
pop_vcf = vcf[wild_dom]
del vcf

# generate ms output
output = open(outputfile_path,'w')

chr_list_vcf = [x[1] for x in pop_vcf.groupby('#CHROM')] # split dataframe into sub dataframe by chr
nLocus = 0
for chr_vcf in chr_list_vcf :
    pos = 0
    max_length = chr_vcf['POS'].max()
    chr_name = chr_vcf['#CHROM'].to_list()[0]
    while pos < max_length:
        sub_chr = chr_vcf[(chr_vcf['POS'] >= pos) & (chr_vcf['POS'] < (pos + locusLength))]
        # convert and write ms file
        head = '//' + chr_name + ':' + str(pos) + '-' + str(pos + locusLength) + '\n'
        output.write(head)
        output.write('segsites: {nseg}\n'.format(nseg=sub_chr.shape[0]))
        segsite = ' '.join([str(round((1.0*i-pos)/locusLength,4)) for i in sub_chr['POS']])
        output.write('positions: {segsite}\n'.format(segsite=segsite))
        for ind in popA.to_list() + popB.to_list():
            line = ''.join([str(x) for x  in sub_chr[ind]])
            output.write(line + '\n')
        output.write('\n')

        nLocus += 1
        pos+=locusLength 

# bpfile
bpfile = open(bpfile_path,'w')
bpfile_data = []
bpfile_data.append('spA={spA} spB={spB} Nref={Nref} mu={mu}\n'.format(spA=nameA,spB=nameB,Nref=Nref,mu=mu))
bpfile_data.append('\t'.join([ str(locusLength) for x in range(nLocus)])+'\n')
bpfile_data.append('\t'.join([ str(len(popA)) for x in range(nLocus)])+'\n')
bpfile_data.append('\t'.join([ str(len(popB)) for x in range(nLocus)])+'\n')
bpfile_data.append('\t'.join([ str(mu*4*Nref*locusLength)  for x in range(nLocus)])+'\n')
bpfile_data.append('\t'.join([ str(rho_over_theta *mu*4*Nref*locusLength)  for x in range(nLocus)])+'\n')
for line in bpfile_data:
    bpfile.write(line)
bpfile.close()

    
