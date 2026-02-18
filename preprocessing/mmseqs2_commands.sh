#!/usr/bin/env bash
# MMseqs2 redundancy removal (run separately for each dataset)
# Thresholds required by the course: 30% identity and 40% coverage

mmseqs createdb INPUT.fasta DB
mmseqs cluster DB CLU TMP --min-seq-id 0.3 -c 0.4 --cov-mode 0
mmseqs createseqfiledb DB CLU REPDB
mmseqs result2flat DB DB REPDB OUTPUT_rep.fasta --use-fasta-header

