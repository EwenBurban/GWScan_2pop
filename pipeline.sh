#!/bin/bash
#SBATCH --job-name=pipe1
#SBATCH -n 30
#SBATCH --mem 10G

### default values #####
DIRSCRIPT=~/DomIsol/script/
PHASELIST="Format Filtering Phasing haplotyping vcftools_analyse pixy_analyse convert2fasta"
PHASES=( $PHASELIST )
PHASE_TO_DO=${PHASES[0]}
FORMAT=linear
C=30
#### functions customs ##
function next_phase() {
	# $1 blockname
	# $2 format of progression
	# $... all phases 
	if [[ $2 == "seq" ]]; then
		return 0
	elif [[ $2 == "linear" ]]; then
		I=2
		arg=( $@ )
		while [[ $1 != ${arg[$I]} ]]; do 
			I=$(($I+1))
		done
		echo ${arg[$(($I+1))]}
	fi
}

### options ######
optstring=":hfd:s:c:"
while getopts ${optstring} arg; do
	case ${arg} in 
		h)
			echo help in coming in further version
			;;
		d)
			DIR="${OPTARG}"
			;;
		s)
			PHASE_TO_DO=${OPTARG}
			;;
		c)
			C=${OPTARG}
			;;
		f)
			FORMAT=seq
			;;
		?)
			echo "Invalid option: -${OPTARG}."
			exit 0
			;;
	esac
done

###################
# DIR is the directory where data are data are strucured like this 
# rawdata --| vcf_files --| all vcf files
#	  		| wildpop.txt
#	  		| dompop.txt
#			| pop.txt 		NB : the popfile contain sample name and the pop at which they come from
#	  		| genome_ref.fasta
#	  		| rho_maps --| all rho map.txt
#
# processed_data --|
# supp_data --|
# .cache

echo the pointed directory is : $DIR
echo the number of core to use is : $C
. /local/env/envconda.sh

#. /local/env/envr-3.6.2.sh
#rm slurm*
rm "$DIR".cache/*
ls "$DIR"rawdata/vcf_files | cat > "$DIR".cache/list_vcf_files.txt
echo $(ls "$DIR".cache)
BLOCKNAME="Format"  
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	# format raw data to keep only GT data
	conda activate bcftools 
	. /local/env/envparallel-20190122.sh
	ls "$DIR"rawdata/vcf_files | cat > "$DIR".cache/list_vcf_files.txt
	mkdir -p "$DIR"processed_data/format_vcf_files
	while IFS= read -r line; do
		echo "bcftools annotate -x 'FORMAT' "$DIR"rawdata/vcf_files/"$line" | cat > "$DIR"processed_data/format_vcf_files/"$line" " >>  "$DIR".cache/format_exec.sh
	done < "$DIR".cache/list_vcf_files.txt
	parallel -a "$DIR".cache/format_exec.sh -j $C
	conda deactivate
	
	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi

BLOCKNAME="Filtering"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	# filter data to keep only bi-allelic sites
	conda activate bcftools
	. /local/env/envparallel-20190122.sh
	mkdir -p "$DIR"processed_data/filtered_vcf_files
	while IFS= read -r line; do
		echo "bcftools view --max-alleles 2 --exclude-types indels "$DIR"processed_data/format_vcf_files/"$line" | cat > "$DIR"processed_data/filtered_vcf_files/"$line" " >> "$DIR".cache/filter_exec.sh
	done < "$DIR".cache/list_vcf_files.txt
	parallel -a "$DIR".cache/filter_exec.sh -j $C
	conda deactivate
	 
	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi

BLOCKNAME="Phasing"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	### phasage
	conda activate shapeit
	mkdir -p "$DIR"processed_data/phased_vcf_files
	while IFS= read -r line; do
		#reg=$(echo $(tail -n1 "$DIR"processed_data/filtered_vcf_files/"$line") | cut -d ' ' -f1 )
		shapeit --input-vcf "$DIR"processed_data/filtered_vcf_files/"$line" -O "$DIR"processed_data/phased_vcf_files/"$line" -T $C 
		shapeit -convert --input-haps "$DIR"processed_data/phased_vcf_files/"$line" --output-vcf "$DIR"processed_data/phased_vcf_files/"$line"
		#shapeit4 --input "$DIR"processed_data/filtered_vcf_files/"$line" --output "$DIR"processed_data/phased_vcf_files/"$line" --thread $C --region $reg
	done < "$DIR".cache/list_vcf_files.txt
	conda deactivate
	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi

#tar -cvzf "$DIR"processed_data/filtered_vcf_files.tar.gz "$DIR"processed_data/filtered_vcf_files/
BLOCKNAME="haplotyping"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	### haplotyping phased data 
	. /local/env/envparallel-20190122.sh
	. /local/env/envperl-5.26.2.sh
	mkdir -p "$DIR"processed_data/haplotype_vcf_files
	while IFS= read -r line; do
		echo "$DIRSCRIPT"make_haplo_vcf.pl "$DIR"processed_data/phased_vcf_files/"$line" "$DIR"processed_data/haplotype_vcf_files/"$line" 1 >> "$DIR".cache/haplotyping_exec.sh
	done < "$DIR".cache/list_vcf_files.txt
	parallel -a "$DIR".cache/haplotyping_exec.sh -j $C
	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi

BLOCKNAME="vcftools_analyse"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	### classic demographic (vcftools)
	conda activate vcftools
	. /local/env/envparallel-20190122.sh
	
	while IFS= read -r line; do
		echo vcftools --vcf "$DIR"processed_data/haplotype_vcf_files/"$line" --weir-fst-pop "$DIR"rawdata/wildpop.txt --weir-fst-pop "$DIR"rawdata/dompop.txt --out "$DIR"processed_data/"$line" >> "$DIR".cache/classic_demo_exec.sh
		echo vcftools --vcf "$DIR"processed_data/haplotype_vcf_files/"$line" --weir-fst-pop "$DIR"rawdata/wildpop.txt --weir-fst-pop "$DIR"rawdata/dompop.txt  --fst-window-size 100000 --out "$DIR"processed_data/"$line" >> "$DIR".cache/classic_demo_exec.sh
		echo vcftools --vcf "$DIR"processed_data/haplotype_vcf_files/"$line" --SNPdensity 100000 --out "$DIR"processed_data/"$line" >> "$DIR".cache/classic_demo_exec.sh
		echo vcftools --vcf "$DIR"processed_data/haplotype_vcf_files/"$line" --TajimaD 100000 --keep "$DIR"rawdata/wildpop.txt --out "$DIR"processed_data/"$line"_wild >> "$DIR".cache/classic_demo_exec.sh
		echo vcftools --vcf "$DIR"processed_data/haplotype_vcf_files/"$line" --TajimaD 100000 --keep "$DIR"rawdata/dompop.txt --out "$DIR"processed_data/"$line"_dom >> "$DIR".cache/classic_demo_exec.sh
		echo vcftools --vcf "$DIR"processed_data/haplotype_vcf_files/"$line" --site-pi --window-pi 100000 --keep "$DIR"rawdata/wildpop.txt --out "$DIR"processed_data/"$line"_wild >> "$DIR".cache/classic_demo_exec.sh
		echo vcftools --vcf "$DIR"processed_data/haplotype_vcf_files/"$line" --site-pi --window-pi 100000 --keep "$DIR"rawdata/dompop.txt --out "$DIR"processed_data/"$line"_dom >> "$DIR".cache/classic_demo_exec.sh
	done < "$DIR".cache/list_vcf_files.txt
	parallel -a "$DIR".cache/classic_demo_exec.sh -j $C
	conda deactivate
	
	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi

BLOCKNAME="pixy_analyse"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	conda activate pixy
	. /local/env/envparallel-20190122.sh
	mkdir -p "$DIR"processed_data/pixy_format
	while IFS= read -r line; do
		echo "bgzip -c "$DIR"processed_data/haplotype_vcf_files/"$line" | cat > "$DIR"processed_data/pixy_format/"$line".gz" >> "$DIR".cache/gzip_exec.sh
		echo tabix "$DIR"processed_data/pixy_format/"$line".gz >> "$DIR".cache/tabix_exec.sh
		echo pixy --stats pi fst dxy --vcf "$DIR"processed_data/pixy_format/"$line".gz --populations "$DIR"rawdata/pop.txt --window_size 100000 --output_prefix "$line"  --bypass_invariant_check 'yes' >> "$DIR".cache/pixy_exec.sh
	done < "$DIR".cache/list_vcf_files.txt
	parallel -a "$DIR".cache/gzip_exec.sh -j $C
	parallel -a "$DIR".cache/tabix_exec.sh -j $C
	parallel -a "$DIR".cache/pixy_exec.sh -j $C
	conda deactivate
	
	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi

BLOCKNAME="convert2fasta"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	# The first objectif is too split the genome ref into locus ( independent locus or just locus )
	conda activate bcftools
	. /local/env/envparallel-20190122.sh
	. /local/env/envperl-5.26.2.sh
	mkdir -p "$DIR"processed_data/sample_fasta
	"$DIRSCRIPT"transform_genome.pl "$DIR"rawdata/genome_ref.fasta "$DIR"processed_data/transformed_genome_ref.fasta
	"$DIRSCRIPT"compil_vcf.pl "$DIR"processed_data/haplotype_vcf_files/ "$DIR"processed_data/all_chr_haplo.vcf
	bgzip -c "$DIR"processed_data/all_chr_haplo.vcf | cat > "$DIR"processed_data/all_chr_haplo.vcf.gz
#	rm "$DIR"processed_data/all_chr_haplo.vcf
	bcftools index "$DIR"processed_data/all_chr_haplo.vcf.gz
	"$DIRSCRIPT"split_genome2locus.pl "$DIR"processed_data/transformed_genome_ref.fasta "$DIR"processed_data/ref_locus.fasta
	while IFS= read -r line; do
		echo bcftools consensus -f "$DIR"processed_data/ref_locus.fasta -s ${line::-1}  "$DIR"processed_data/all_chr_haplo.vcf.gz -o "$DIR"processed_data/sample_fasta/"${line::-1}"_dom.fasta >> "$DIR".cache/consensus_exec.sh
		echo "$DIRSCRIPT"header_edit.pl  "$DIR"processed_data/sample_fasta/"${line::-1}"_dom.fasta "$DIR"processed_data/samples.fasta dom ${line::-1} >> "$DIR".cache/header_edit_exec.sh
	done < "$DIR"rawdata/dompop.txt

	while IFS= read -r line; do
		echo bcftools consensus -f "$DIR"processed_data/ref_locus.fasta -s ${line::-1}  "$DIR"processed_data/all_chr_haplo.vcf.gz -o "$DIR"processed_data/sample_fasta/"${line::-1}"_wild.fasta >> "$DIR".cache/consensus_exec.sh
		echo "$DIRSCRIPT"header_edit.pl  "$DIR"processed_data/sample_fasta/"${line::-1}"_wild.fasta "$DIR"processed_data/samples.fasta wild ${line::-1} >> "$DIR".cache/header_edit_exec.sh
	done < "$DIR"rawdata/wildpop.txt
	parallel -a "$DIR".cache/consensus_exec.sh -j $C
	parallel -a "$DIR".cache/header_edit_exec.sh -j 1
	conda deactivate
	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi

BLOCKNAME="oldblock"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	# The first objectif is too split the genome ref into locus ( independent locus or just locus )
	conda activate bcftools
	. /local/env/envparallel-20190122.sh
	. /local/env/envperl-5.26.2.sh
	mkdir -p "$DIR"processed_data/consensus_fasta
	"$DIRSCRIPT"transform_genome.pl "$DIR"rawdata/genome_ref.fasta "$DIR"processed_data/transformed_genome_ref.fasta
	"$DIRSCRIPT"compil_vcf.pl "$DIR"processed_data/haplotype_vcf_files/ "$DIR"processed_data/all_chr_haplo.vcf
	bgzip -c "$DIR"processed_data/all_chr_haplo.vcf | cat > "$DIR"processed_data/all_chr_haplo.vcf.gz
#	rm "$DIR"processed_data/all_chr_haplo.vcf
	bcftools index "$DIR"processed_data/all_chr_haplo.vcf.gz
	"$DIRSCRIPT"fullgenome2fasta.pl "$DIR"processed_data/transformed_genome_ref.fasta "$DIR"processed_data/splitted_ref_genome.fasta
	while IFS= read -r line; do
		echo bcftools consensus -f "$DIR"processed_data/splitted_ref_genome.fasta -s ${line::-1}  "$DIR"processed_data/all_chr_haplo.vcf.gz -o "$DIR"processed_data/consensus_fasta/"${line::-1}"_dom.fasta >> "$DIR".cache/consensus_exec.sh
		echo "$DIRSCRIPT"header_edit.pl  "$DIR"processed_data/consensus_fasta/"${line::-1}"_dom.fasta "$DIR"processed_data/consensus.fasta dom ${line::-1} >> "$DIR".cache/header_edit_exec.sh
	done < "$DIR"rawdata/dompop.txt

	while IFS= read -r line; do
		echo bcftools consensus -f "$DIR"processed_data/splitted_ref_genome.fasta -s ${line::-1}  "$DIR"processed_data/all_chr_haplo.vcf.gz -o "$DIR"processed_data/consensus_fasta/"${line::-1}"_wild.fasta >> "$DIR".cache/consensus_exec.sh
		echo "$DIRSCRIPT"header_edit.pl  "$DIR"processed_data/consensus_fasta/"${line::-1}"_wild.fasta "$DIR"processed_data/consensus.fasta wild ${line::-1} >> "$DIR".cache/header_edit_exec.sh
	done < "$DIR"rawdata/wildpop.txt
	parallel -a "$DIR".cache/consensus_exec.sh -j $C
	parallel -a "$DIR".cache/header_edit_exec.sh -j 1
	conda deactivate

	sbatch ~/mDILS/bin/DILS_2pop.sh "$DIR"DILS/genome_scan.yaml

	PHASE_TO_DO=$(next_phase $BLOCKNAME $FORMAT ${PHASES[@]})
fi


BLOCKNAME="globalstat"
if [[ $BLOCKNAME == $PHASE_TO_DO ]]; then
	conda activate bcftools
	. /local/env/envparallel-20190122.sh
	. /local/env/envperl-5.26.2.sh
	mkdir -p "$DIR"processed_data/consensus_fasta
	mkdir -p "$DIR"processed_data/ref_chr
	"$DIRSCRIPT"transform_genome.pl "$DIR"rawdata/genome_ref.fasta "$DIR"processed_data/transformed_genome_ref.fasta
	"$DIRSCRIPT"fullgenome2fasta.pl "$DIR"processed_data/transformed_genome_ref.fasta "$DIR"processed_data/ref_chr/

	if [[ ! -f  "$DIR"processed_data/all_chr_haplo.vcf.gz ]];then
		"$DIRSCRIPT"compil_vcf.pl "$DIR"processed_data/haplotype_vcf_files/ "$DIR"processed_data/all_chr_haplo.vcf
		bgzip -c "$DIR"processed_data/all_chr_haplo.vcf | cat > "$DIR"processed_data/all_chr_haplo.vcf.gz
		bcftools index "$DIR"processed_data/all_chr_haplo.vcf.gz
	fi
	
	ls "$DIR"processed_data/ref_chr  | cat > "$DIR".cache/list_chr_fasta_files.txt
	while IFS= read -r chr;do
		rm  "$DIR"processed_data/consensus_fasta/"${chr::-10}"_allind.fasta
		while IFS= read -r line; do
			echo bcftools consensus -f "$DIR"processed_data/ref_chr/"$chr" -s ${line::-1}  "$DIR"processed_data/all_chr_haplo.vcf.gz -o "$DIR"processed_data/consensus_fasta/"${chr::-10}"_"${line::-1}"_dom.fasta >> "$DIR".cache/consensus_exec.sh
			echo "$DIRSCRIPT"header_edit.pl  "$DIR"processed_data/consensus_fasta/"${chr::-10}"_"${line::-1}"_dom.fasta "$DIR"processed_data/consensus_fasta/"${chr::-10}"_allind.fasta dom ${line::-1} >> "$DIR".cache/header_edit_exec.sh
		done < "$DIR"rawdata/dompop.txt

		while IFS= read -r line; do
			echo bcftools consensus -f "$DIR"processed_data/ref_chr/"$chr" -s ${line::-1}  "$DIR"processed_data/all_chr_haplo.vcf.gz -o "$DIR"processed_data/consensus_fasta/"${chr::-10}"_"${line::-1}"_wild.fasta >> "$DIR".cache/consensus_exec.sh
			echo "$DIRSCRIPT"header_edit.pl  "$DIR"processed_data/consensus_fasta/"${chr::-10}"_"${line::-1}"_wild.fasta "$DIR"processed_data/consensus_fasta/"${chr::-10}"_allind.fasta wild ${line::-1} >> "$DIR".cache/header_edit_exec.sh
		done < "$DIR"rawdata/wildpop.txt

		cp ~/mDILS/genome_scan_config.yaml "$DIR".cache/"${chr::-10}"_scan_mDILSconfig.yaml
		sed -i "s\infile:.*$\infile: ${DIR}processed_data/consensus_fasta/${chr::-10}_allind.fasta\g" "$DIR".cache/"${chr::-10}"_scan_mDILSconfig.yaml 
		sed -i "s\timeStamp:.*$\timeStamp: ${DIR}DILS/${chr::-10}_scan\g" "$DIR".cache/"${chr::-10}"_scan_mDILSconfig.yaml 
		sed -i "s\config_yaml:.*$\config_yaml: ${DIR}.cache/${chr::-10}_scan_mDILSconfig.yaml\g" "$DIR".cache/"${chr::-10}"_scan_mDILSconfig.yaml 

		echo sbatch ~/mDILS/bin/DILS_2pop.sh "$DIR".cache/"${chr::-10}"_scan_mDILSconfig.yaml >> "$DIR".cache/mDILS_exec.sh
		# add gestion des config de mDILS + un exec pour mDILS
	done < "$DIR".cache/list_chr_fasta_files.txt
	parallel -a "$DIR".cache/consensus_exec.sh -j $C
	parallel -a "$DIR".cache/header_edit_exec.sh -j 1
	parallel -a "$DIR".cache/mDILS_exec.sh -j $C
	conda deactivate
fi
