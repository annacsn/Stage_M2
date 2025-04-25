#! /bin/bash

# Script written by Romain Coppee / Anna Cosson
# Creation date: 09/02/2020
# Last modification: 18/02/2025

####--------General Goal: Produce VCF files containing only high-quality SNPs from whole genome sequencing data

# Set up the conda environment for GATK
__conda_setup="$('/opt/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
conda activate /opt/anaconda3/envs/gatk_env/

# Path to the PICARD software
PICARD=/Users/adm-loc/Downloads/picard.jar

# Define variables for BAM files
FILES_BAM=*.bam

# Define variables for BAM index files
FILES_BAI=*.bai

# Define variables for fixed BAM files
FILES_FIX=*.fix

# Path to the reference genome
REF_GEN=/Volumes/ANNA_DD/Projet_Asie_du_Sud_Est/Ref_genome/ref/Pfalciparum.genome.fasta

# Path to known variant sites (for base recalibration)
GEN_CROSS=/Volumes/ANNA_DD/Projet_Asie_du_Sud_Est/Ref_genome/known_sites

# Text file listing the chromosomes to process
CHROM_LIST="chromosomes.txt"

###############################################################################################################
#####-------Goal 1: Complete header information, fix BAM files, and index them

# Step 1: Mark duplicates with PicardTools
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
    echo "✅ Mark Duplicates PROCESSED for $f"
done

# Step 2: Remove original .bam files and rename the fixed ones
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "✅ Renamed fixed BAM after MarkDuplicates: $f"
done

# Step 3: Add or Replace Read Groups with Picard
for f in `ls $FILES_BAM`
do
    file_name="${f##*/}"
    onlyFileName="${file_name%.*}"  # Remove .bam extension for SM tag
    onlyFileName="${onlyFileName%.*}"  # Also remove .sorted if present

    java -jar $PICARD AddOrReplaceReadGroups \
        -I $f \
        -O $f.fix \
        -LB WGA \
        -PL illumina \
        -PU NA \
        -SM $onlyFileName
    echo "✅ AddOrReplaceReadGroups PROCESSED for $f"
done

# Step 4: Again, remove original .bam files and rename the new fixed ones
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "✅ Renamed fixed BAM after AddOrReplaceReadGroups: $f"
done

# Step 5: Index the BAM files with samtools
for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "✅ BAM indexing PROCESSED for $f"
done

# Step 6: Base Quality Score Recalibration (BQSR) - create recalibration tables
for f in `ls $FILES_BAM`
do
    gatk BaseRecalibrator \
        -R $REF_GEN \
        -I $f \
        --known-sites $GEN_CROSS/3d7_hb3.combined.final.vcf.gz \
        --known-sites $GEN_CROSS/hb3_dd2.combined.final.vcf.gz \
        --known-sites $GEN_CROSS/7g8_gb4.combined.final.vcf.gz \
        -O $f.pass1.table 
    echo "✅ BaseRecalibrator PROCESSED for $f"
done

# Step 7: Apply BQSR using the recalibration tables
for f in `ls $FILES_BAM`
do
    gatk ApplyBQSR \
        -R $REF_GEN \
        -I $f \
        -bqsr $f.pass1.table \
        -O $f.fix
done

# Step 8: Remove original BAM files and rename the recalibrated ones
rm $FILES_BAM
for f in `ls $FILES_FIX`
do 
    mv -- "$f" "${f%.fix}"
    echo "✅ Renamed BAM after ApplyBQSR: $f"
done

# Step 9: Re-index the recalibrated BAM files
rm $FILES_BAI
for f in `ls $FILES_BAM`
do
    samtools index $f
    echo "✅ Re-indexed recalibrated BAM: $f"
done

###############################################################################################################
#####-------Goal 2: Perform variant calling (GVCF mode) with GATK HaplotypeCaller

# Step 10: Variant Calling with GATK HaplotypeCaller (per chromosome, parallelized)
for f in $FILES_BAM; do
    sample_name="${f%%.*}" # Extract sample name
    mkdir -p HaplotypeCaller_temp/${sample_name}

    # Launch HaplotypeCaller separately for each chromosome
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
    echo "✅ HaplotypeCaller PARALLEL COMPLETED for $sample_name"

    # Step 11: Combine GVCFs from all chromosomes into a single GVCF per sample

    # Collect all chromosome GVCF files
    gvcf_files=()
    while IFS= read -r file; do
        gvcf_files+=("$file")
    done < <(find HaplotypeCaller_temp/${sample_name} -name "${sample_name}.*.g.vcf" | sort)

    # Check that GVCF files were found
    if [ ${#gvcf_files[@]} -eq 0 ]; then
        echo "❌ No GVCF files found for $sample_name"
        exit 1
    fi

    # (Optional) Display found files
    echo "GVCF files found for $sample_name:"
    for f in "${gvcf_files[@]}"; do echo "  $f"; done

    # Combine GVCFs
    gatk CombineGVCFs \
        -R $REF_GEN \
        $(for gvcf in "${gvcf_files[@]}"; do echo "-V $gvcf"; done) \
        -O ${sample_name}.g.vcf

    echo "✅ CombineGVCFs COMPLETED for $sample_name"
done

# Deactivate conda environment
conda deactivate
