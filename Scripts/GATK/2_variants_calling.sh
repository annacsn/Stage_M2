#! /bin/bash

#Script Romain Coppee & Anna Cosson
#Creation data: 09/02/2020
#Last modification: 18/04/2025

####--------General Goal: Produce VCF files containing only high-quality SNPs from whole genome sequencing data

###############################################################################################################
#####-------Preparation/location of softwares and materials
#Adding samtools/bcftools to the PATH environment
#export PATH="$PATH:/usr/bin/bcftools-1.9"
#export PATH="$PATH:/usr/bin/samtools-1.9"
#export PATH="$PATH:/usr/bin/htslib-1.9"
#source ~/.profile

#Adding GATK to the PATH environment
#export PATH=/home/virologie/Documents/gatk-4.1.8.1/:$PATH

__conda_setup="$('/opt/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate /opt/anaconda3/envs/gatk_env/

#location of PICARD software
PICARD=/Users/adm-loc/Downloads/picard.jar

#FILES contains the BAM file for each sample
FILES_BAM=*.bam

#FILES contains the indexed BAM file for each sample
FILES_BAI=*.bai

#FILES contains the fixed, indexed BAM file for each sample
FILES_FIX=*.fix

#Location of the Reference genome
REF_GEN=/Volumes/ANNA_DD/Projet_Asie_du_Sud_Est/Ref_genome/ref/Pfalciparum.genome.fasta

#Location of genetic crosses
GEN_CROSS=/Volumes/ANNA_DD/Projet_Asie_du_Sud_Est/Ref_genome/known_sites

CHROM_LIST="chromosomes.txt"


###############################################################################################################
#####-------Goal 1: Complete Header to fix BAM files, then index them



#Mark Duplicates with PicardTools
for f in `ls $FILES_BAM`
do
    java -jar $PICARD MarkDuplicates \
        -INPUT $f \
        -REMOVE_DUPLICATES true \
        -VALIDATION_STRINGENCY LENIENT \
        -AS true \
        -METRICS_FILE metrics \
        -CREATE_INDEX true \
        -OUTPUT $f.fix
    echo "Mark Duplicates PROCESSED $f"
done

#first remove all .bam files, and rename the fixed .bam.fix to .bam
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing mark duplicates PROCESSED $f"
done



#complete the header of each bam file using the AddOrReplaceReadGroups from Picard 
for f in `ls $FILES_BAM`
do
    file_name="${f##*/}"
    onlyFileName="${file_name%.*}" #remove.bam for SM information
    onlyFileName="${onlyFileName%.*}" #remove .sorted SM information
    java -jar $PICARD AddOrReplaceReadGroups \
        -I $f \
        -O $f.fix \
        -LB WGA \
        -PL illumina \
        -PU NA \
        -SM $onlyFileName
    echo "AddOrReplaceReadGroups PROCESSED $f"
done


#first remove all .bam files, and rename the fixed .bam.fix to .bam
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing AddOrReplaceReadGroups PROCESSED $f"
done

#index each fixed bam file using samtools
for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "index PROCESSED $f"
done


#Recalibrate base qualities in a .bam file so that quality metrics match actual observed error rates
for f in `ls $FILES_BAM`
do
    gatk BaseRecalibrator \
        -R $REF_GEN \
        -I $f \
        --known-sites $GEN_CROSS/3d7_hb3.combined.final.vcf.gz \
        --known-sites $GEN_CROSS/hb3_dd2.combined.final.vcf.gz \
        --known-sites $GEN_CROSS/7g8_gb4.combined.final.vcf.gz \
        -O $f.pass1.table 
    echo "BaseRecalibrator PROCESSED $f"
done

#Apply BQSR to input BAM file
for f in `ls $FILES_BAM`
do
    gatk ApplyBQSR \
        -R $REF_GEN \
        -I $f \
        -bqsr $f.pass1.table \
        -O $f.fix
done

#first remove all .bam files, and rename the fixed .bam.fix to .bam
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "fixing PrintReads PROCESSED $f"
done

#remove the previous indexed bam files and create new indexed files
rm $FILES_BAI
for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "index PROCESSED $f"
done

#generate a complete variant calling from each fixed bam file
# Étape 1 : Appel parallèle HaplotypeCaller pour chaque chromosome
for f in $FILES_BAM; do
    sample_name="${f%%.*}"
    mkdir -p HaplotypeCaller_temp/${sample_name}

    while read chr; do
        gatk HaplotypeCaller \
            -R $REF_GEN \
            -I $f \
            -ERC GVCF \
            --native-pair-hmm-threads 2 \
            --max-alternate-alleles 6 \
            -L $chr \
            -O HaplotypeCaller_temp/${sample_name}/${sample_name}.${chr}.g.vcf &
    done < $CHROM_LIST

    wait
    echo "HaplotypeCaller PARALLEL COMPLETED for $sample_name"

    # Étape 2 : Combine les fichiers GVCF (1 par chromosome) en un seul par échantillon

    # Liste les fichiers g.vcf
    gvcf_files=()
    while IFS= read -r file; do
        gvcf_files+=("$file")
    done < <(find HaplotypeCaller_temp/${sample_name} -name "${sample_name}.*.g.vcf" | sort)

    # Vérifier s'il y a bien des fichiers à combiner
    if [ ${#gvcf_files[@]} -eq 0 ]; then
        echo "❌ Aucun fichier GVCF trouvé pour $sample_name"
        exit 1
    fi

    # (optionnel) Affiche les fichiers trouvés
    echo "Fichiers GVCF trouvés pour $sample_name :"
    for f in "${gvcf_files[@]}"; do echo "  $f"; done

    # Combine les fichiers
    gatk CombineGVCFs \
        -R $REF_GEN \
        $(for gvcf in "${gvcf_files[@]}"; do echo "-V $gvcf"; done) \
        -O ${sample_name}.g.vcf

echo "✅ CombineGVCFs terminé pour $sample_name"



done

conda deactivate
