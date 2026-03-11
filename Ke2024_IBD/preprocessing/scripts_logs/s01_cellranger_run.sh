#!/bin/bash

### config
FASTQ_DIR="../../data/fastq" # dir with fastq subdirs, sample names as subdir names
TRANSCRIPTOME="/home/cellranger/refdata-gex-GRCh38-2024-A"
OUT_DIR="../../data/cellranger"
NUM_THREADS=30
LOGFILE="cellranger_run_$(date '+%Y%m%d_%H%M%S').log"

### logging
echo "🚀 Starting Cell Ranger at $(date)" | tee -a "$LOGFILE"
mkdir -p "$OUT_DIR"

### loop through dirs
for SAMPLE_PATH in "$FASTQ_DIR"/*/; do
    SAMPLE_NAME=$(basename "$SAMPLE_PATH")
    echo "🤖 Processing sample: $SAMPLE_NAME at $(date)" | tee -a "$LOGFILE"
    echo "📂 FASTQ path: $SAMPLE_PATH" | tee -a "$LOGFILE"

    START_TIME=$(date +%s)

    cellranger count \
        --id="$SAMPLE_NAME" \
        --transcriptome="$TRANSCRIPTOME" \
        --fastqs="$SAMPLE_PATH" \
        --sample="$SAMPLE_NAME" \
        --localcores="$NUM_THREADS" \
        --create-bam=true \
        --include-introns=true \
        --localmem=$((NUM_THREADS * 7)) \
        --output-dir="$OUT_DIR/$SAMPLE_NAME" \
        &>> "$LOGFILE"

    END_TIME=$(date +%s)
    DURATION=$((END_TIME - START_TIME))
    echo "✅ Finished $SAMPLE_NAME at $(date) ⏳ Duration: ${DURATION}s" | tee -a "$LOGFILE"
    echo "❤️🧡💛💚💙💜❤️🧡💛💚💙💜❤️" | tee -a "$LOGFILE"
done

### wohoo
echo "🎉 All finished at $(date)" | tee -a "$LOGFILE"
echo "💾 Log saved to $LOGFILE"