#!/bin/bash
for folder in KCIS*;do
    echo "Sample: $folder"
    zcat "$folder/matrix.mtx.gz" | head -n 10 | cut -d ' ' -f 1-4
    echo ""
done
