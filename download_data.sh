#!/bin/bash
set -euo pipefail

# 1. Create a local data folder on your host machine
#mkdir -p data
cd data

echo "[INFO] Commencing file downloads to local host..."

# 2. Download the CRAM files and reference genome locally
wget -c https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/COLO829T_TEST.cram
wget -c https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/COLO829T_TEST.cram.crai
wget -c https://aveit.s3.us-east-1.amazonaws.com/misc/INTERVIEW/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

cd ..

echo "[INFO] Running samtools faidx inside the container to create the .fai file..."

# 3. Mount your local ./data folder to the container's /workspace/data/ path
# This maps the files perfectly and generates the .fai file right next to your FASTA
docker run --rm \
  -v "$(pwd)/data:/workspace/data" \
  cram-pipeline:latest \
  samtools faidx /workspace/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa

echo "[SUCCESS] All files are now available under /workspace/data/ during runtime!"

