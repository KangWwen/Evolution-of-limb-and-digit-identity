#!/bin/bash
GENOME_FA="~/ref/Gallus_gallus.GRCg6a.dna.toplevel.fa"
GTF_FILE="~/ref/Gallus_gallus.Gallus_gallus-5.0.91.gtf"
INDEX_PREFIX="~/ref/chick_mRNA_index/gallus_tran"

READS_DIR="./Reads"
BAM_DIR="./BamFiles"
COUNT_DIR="./CountFiles"
LOG_DIR="./Logs"

mkdir -p ${BAM_DIR}
mkdir -p ${COUNT_DIR}
mkdir -p ${LOG_DIR}

THREADS=2


if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo "Building HISAT2 index..."
    hisat2-build -p ${THREADS} ${GENOME_FA} ${INDEX_PREFIX} \
        2>&1 | tee ${LOG_DIR}/hisat2_build.log
    echo "Index building completed!"
else
    echo "Index already exists, skipping..."
fi

for READ1 in ${READS_DIR}/*_paired_1.fq.gz; do
    if [ ! -f "${READ1}" ]; then
        echo "No read files found in ${READS_DIR}"
        exit 1
    fi
    BASENAME=$(basename ${READ1} _paired_1.fq.gz)
    READ2="${READS_DIR}/${BASENAME}_paired_2.fq.gz"

    if [ ! -f "${READ2}" ]; then
        echo "Warning: Paired file ${READ2} not found, skipping ${BASENAME}"
        continue
    fi
    
    SAM_FILE="${BAM_DIR}/${BASENAME}.sam"
    BAM_FILE="${BAM_DIR}/${BASENAME}.bam"
    SORTED_BAM="${BAM_DIR}/${BASENAME}_sorted.bam"

    echo "Running HISAT2 alignment and converting to BAM..."
    hisat2 -p ${THREADS} \
        --dta \
        -x ${INDEX_PREFIX} \
        -1 ${READ1} \
        -2 ${READ2} \
        2> ${LOG_DIR}/${BASENAME}_hisat2.log | \
        samtools view -@ ${THREADS} -bS - | \
        samtools sort -@ ${THREADS} -o ${SORTED_BAM} -
    
    # Index the BAM file
    echo "Indexing BAM file..."
    samtools index ${SORTED_BAM}
    
    echo "Alignment completed for ${BASENAME}"
done

for BAM_FILE in ${BAM_DIR}/*_sorted.bam; do
    if [ ! -f "${BAM_FILE}" ]; then
        echo "No BAM files found in ${BAM_DIR}"
        exit 1
    fi
    
    BASENAME=$(basename ${BAM_FILE} _sorted.bam)
    COUNT_FILE="${COUNT_DIR}/${BASENAME}_count.txt"
    
    echo ""
    echo "Counting features for: ${BASENAME}"
    
    featureCounts \
        -T ${THREADS} \
        -p \
        -t exon \
        -g gene_id \
        -a ${GTF_FILE} \
        -o ${COUNT_FILE} \
        ${BAM_FILE} \
        2>&1 | tee ${LOG_DIR}/${BASENAME}_featureCounts.log
    
    echo "Feature counting completed for ${BASENAME}"
done

echo ""
echo "=========================================="
echo "Pipeline Completed Successfully!"
echo "=========================================="
echo ""
echo "Summary:"
echo "- BAM files saved in: ${BAM_DIR}"
echo "- Count files saved in: ${COUNT_DIR}"
echo "- Log files saved in: ${LOG_DIR}"
echo ""
echo "Total samples processed: $(ls ${COUNT_DIR}/*_count.txt 2>/dev/null | wc -l)"
echo ""