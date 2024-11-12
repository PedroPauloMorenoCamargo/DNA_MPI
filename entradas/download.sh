#!/bin/bash

# Loop over chromosome numbers from 1 to 22
for X in {1..22}
do
    # Construct the URL
    url="ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/snp147Mask/chr${X}.subst.fa.gz"
    
    # Download the file
    echo "Downloading chr${X}.subst.fa.gz..."
    wget "$url"
    
    # Unzip the file
    echo "Unzipping chr${X}.subst.fa.gz..."
    gunzip "chr${X}.subst.fa.gz"
    
    echo "chr${X}.subst.fa downloaded and unzipped."
done

echo "All downloads and unzips complete."
