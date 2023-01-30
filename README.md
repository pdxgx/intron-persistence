# intronomer
Tool for calculating intron persistence, as defined in https://doi.org/10.1186/s13059-022-02789-6

### Usage:
Python 3 code to calculate intron persistence for a long-read RNA-seq sample

1. Assigns best matches between aligned pacbio reads and annotated transcripts
2. Calculates intron persistence values from read-transcript sets

Sample run:
```
python intronomer/intronomer/intronomer.py -g ANNOTATION_FILE -a ALIGNED_READS_FILE -p PROJECT_FLAG -o OUTPUT_DIRECTORY -t PREPROCESSED_READ-TX_FILE
```

#### Required parameters
- ALIGNED_READS_FILE (.bam or .sam format)
    - output of long RNA-seq reads mapped to a reference with a splice-aware
      aligner
    - If using PacBio reads, full length non-concatamer poly-a selected long
      reads are recommended; please do not use reads that have been clustered
      and collapsed into a consensus isoform.
- ANNOTATION_FILE (.gtf format)
    - annotation based on the same reference genome to which the long RNA-seq
      reads were mapped. Any potentially retained persistent introns will be
      defined with respect to this annotation.

#### Optional parameters
- PROJECT_FLAG
    - a flag or project name to add to output filenames
- OUTPUT_DIRECTORY
    - desired output directory; if not given, the alignment file directory
      will be used to write new files
- PREPROCESSED_READ-TX_FILE
    - if the RI_txs_to_read_ids file was previously generated (first step in
      the calculation) you can add it here to avoid re-generating this.

