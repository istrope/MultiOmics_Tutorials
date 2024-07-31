#!/bin/bash
for folder in KCIS*;do
    echo "Sample: $folder"
    zcat "$folder/features.tsv.gz" | head -n 2
    echo ""
done
