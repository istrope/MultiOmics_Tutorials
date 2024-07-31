#!/bin/bash
for folder in KCIS*;do
    echo "Sample: $folder"
    zcat "$folder/barcodes.tsv.gz" | head -n 2
    echo ""
done
