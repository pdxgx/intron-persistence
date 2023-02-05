# intronomer
Tool for calculating intron persistence, as defined in https://doi.org/10.1186/s13059-022-02789-6

### Overview
Python 3 code to calculate intron persistence for a long-read RNA-seq sample

1. Assigns best matches between aligned pacbio reads and annotated transcripts
2. Calculates intron persistence values from read-transcript sets

### Requirements
- `python3` (tested with v3.8.1)
- `pysam` (tested with v0.16.0.1, using samtools v1.10)
- `scipy` (tested with v1.5.4)
- `pandas` (tested with v1.1.4)
- `numpy` (tested with v1.19.4)

### Usage

#### Sample run:
```
python intronomer/intronomer/intronomer.py -g ANNOTATION_FILE -a ALIGNED_READS_FILE -p PROJECT_FLAG -o OUTPUT_DIRECTORY -t PREPROCESSED_READ-TX_FILE
```

#### Required parameters
- `ALIGNED_READS_FILE` (.bam or .sam format)
    - output of long RNA-seq reads mapped to a reference with a splice-aware aligner
    - If using PacBio reads, full length non-concatamer poly-a selected long reads are recommended; please do not use reads that have been clustered and collapsed into a consensus isoform.
- `ANNOTATION_FILE` (.gtf format)
    - annotation based on the same reference genome to which the long RNA-seq reads were mapped. Any potentially retained persistent introns will be defined with respect to this annotation.

#### Optional parameters
- `PROJECT_FLAG`
    - a flag or project name to add to output filenames
- `OUTPUT_DIRECTORY`
    - desired output directory; if not given, the alignment file directory will be used to write new files
- `PREPROCESSED_READ-TX_FILE`
    - if the "RI_txs_to_read_ids" file was previously generated (first step in the calculation) you can add it here to avoid re-generating this.

### Setup

#### Virtual environment

You may run `intronomer` by setting up a virtual environment with the provided .yml file. With conda, make the new environment with:

```
conda env create -f ./yml/intronomer.yml
```

#### Containerization

You may run `intronomer` from a Docker container by setting up a new image using the provided Dockerfile. Do this by navigating to 
the `intronomer` folder and run the following:

```
docker build -t intronomer:latest .
```

We specified the new image name `intronomer:latest` using the tag `-t`. See `docker run --help` for details.

Now, running `docker images` should show something like the following:

```
REPOSITORY              TAG       IMAGE ID       CREATED         SIZE  
intronomer              latest    22b7effb653e   6 minutes ago   580MB
```

Now, instead of navigating to a virtual environment or another directory, `intronomer` can be run from your current directory using:

```
docker run -w $PWD -v $PWD:$PWD intronomer:latest -g ANNOTATION_FILE -a ALIGNED_READS_FILE -p PROJECT_FLAG -o OUTPUT_DIRECTORY -t PREPROCESSED_READ-TX_FILE
```

This used the flags `-w` to specify the container working directory, `-v` to mount the local path so that it is visible to the container, and `--rm` to remove the container on exit, which is good practice. See `docker run --help` for details.

### Runnable example

This runnable example demonstrates how to use `intronomer` and describes its outputs. 

#### Example files

Files for the runnable example have been provided in the `./exe/` directory, and these include:

* `exe.bam/exe.bam.bai` : Small subset of 1000 mapped long reads from the human stem cell sample SRP098984/SAMN07611993.
* `gencode.v35.annotation.subset.gtf` : Small sample annotation from GENCODE v35 (chr###).

#### Compute persistence

Process the example data with `intronomer` using:

```
python ./src/intronomer.py -g './exe/gencode.v35.annotation.exe.gtf' -a './exe/exe.bam' -o 'example' -p 'new-example'
```

If you already set up the Docker image (see above), you can instead navigate to the `/exe/` directory and use following (note, `sudo` may be required):

```
docker run -w $PWD -v $PWD:$PWD --rm intronomer:latest -g 'gencode.v35.annotation.exe.gtf' -a 'exe.bam' -o 'example' -p 'new-example' -rm
```

On run success, you should see a message summarizing the results:

```
7 transcripts with at least 5 reads and some persistent RIs
```

There should also be a new folder called `./example/` containing two results tables, a .tsv called something like `RI_txs_to_read_ids_new-example_.tsv` which contains read-transcript mapping info, and a .csv called something like `new-example_intron-persistence_02-04-2023_14.26.08.csv` containing persistences computed by intron (details below). The datetimes in the new table filenames will vary. Also, note how our project flag was included in each table's filename.

#### Output file descriptions

##### Persistence table

The .csv table called `"intron-persistence_*"` contains the computed persistences organized by intron feature details, including genome coordinates. The column names in this .csv are:

1. `intron` : The intron's uniquely defining name, or in other words its genome coordinates.
2. `transcript` : Transcript ID of transcript containing the intron, denoted as en Ensembl ID with ENST*.
3. `persistence` : Computed persistence of the intron.
4. `position` : Relative transcript position. This is the decimal unit value (e.g. 0-1) of the intron's start location in relation to the transcript starting coordinate.

##### Read-transcript mappings table

The .tsv file called `"RI_txs_to_read_ids*"` contains mappings of reads from the .bam file to transcripts in the .gtf annotation file. The column names of the .tsv file are:

1. `transcript` : Transcript ID of transcript containing the intron, denoted as en Ensembl ID with ENST*.
2. `read_id` : Numeric index of the line for the read in the provided .bam file.
3. `match_type` : Flag denoting status of the read-transcript mapping. Can be either of:
    * `skipped_splicing-partial`
    * `skipped_splicing-full_length`
    * `all_introns-partial`
    * `all_introns-full_length`
4. `region_tx_intron_count` : Total introns in the transcript, from the provided .gtf.
5. `unspliced_intron_count` : Total unspliced introns detected in the transcript.
6. `orig_read_id` : The full-length read ID in the .bam file.
7. `read_length` : Length of the read, in bases.
8. `read_left` : Upstream genome coordinate position of the read.
9. `read_right` : Downstream genome coordinate position of the read.
10. `chrom` : Chromosome for the read and transcript.
11. `strand` : Strand for the read and transcript, either "-" for minus strand or "+" for positive strand.
12. `read_introns` : Full list of the corrdinates of introns fully spanned by the read, separated by semicolons (";").
13. `tx_length` : Transcript length in bases.
14. `tx_left` : Upstream genome coordinate position of the transcript.
15. `tx_right` : Downstream genome coordinate position of the transcript.
16. `tx_introns` : Full list of coordinates for introns covered by the transcript, separated by semicolors (";").
17. `RI-read_ratio`  : Fraction of the total transcript covered by the read.