#!/bin/bash

set -e

EMU_REF="${HOME}/ref/ZJU2"
MOUSE_REF="${HOME}/ref/GRCm38"
DATA_DIR="${HOME}/00.data/forelimbs"
OUTPUT_DIR="./CellRanger_Results"

mkdir -p ${OUTPUT_DIR}

cd ${EMU_REF}
cellranger mkgtf ZJU2.gtf ZJU2.filter.gtf --attribute=gene_biotype:protein_coding
cellranger mkref --genome=ZJU2 --fasta=ZJU2.fa --genes=ZJU2.filter.gtf
mv ZJU2 cellrangerIndex

cd ${MOUSE_REF}
cellranger mkgtf GRCm38.gtf GRCm38.filter.gtf --attribute=gene_biotype:protein_coding
cellranger mkref --genome=GRCm38 --fasta=GRCm38.fa --genes=GRCm38.filter.gtf
mv GRCm38 cellrangerIndex

cd ${OUTPUT_DIR}
cellranger count --id=EF25 --transcriptome=${EMU_REF}/cellrangerIndex --fastqs=${DATA_DIR}/Emu_HH25
cellranger count --id=MF105 --transcriptome=${MOUSE_REF}/cellrangerIndex --fastqs=${DATA_DIR}/Mouse_E105

echo "Completed!"