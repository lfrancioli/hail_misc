#!/bin/bash

python extract_introns.txt.py --input ucsc_gencode_v19.gz --output ucsc_gencode_v19.introns.bed
bedtools getfasta -fi ~/resources/Homo_sapiens_assembly19.fasta -bed ucsc_gencode_v19.introns.bed -name -bedOut > ucsc_gencode_v19.introns.fasta.bed