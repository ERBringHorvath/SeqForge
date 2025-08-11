# **<ins>SeqForge<ins/>**

## Design Philosophy:

SeqForge emphasizes clarity, flexibility, and scale. Each module is standalone but interoperable with others, allowing researchers to plug in just what they need or use the full pipeline. Designed for power users but accessible enough for those just entering the field, SeqForge aims to grow as the field of microbial genomics and bioinformatics expands.

### Note on FASTA file submission:

All FASTA-handling modules support the following inputs:
*    Single FASTA file (.fa, .fna, .faa, .ffn, .fas, .fasta)
*    Single compressed FASTA file (.gz, .zip)
*    Directory of FASTA files
*    Archive of FASTA files (.zip, .tar, .tar.gz, .tgz)

### **SeqForge is currently divided into three modules:**

### <ins>Module 1: Genome Search<ins/>

Purpose: Rapid database creation and high-throughput querying <br/>
*   makedb<br/>
    Create BLAST-compatible databases (makeblastdb). Supports both nucleotide and protein databases, with multiprocessing support to boost performance on large datasets.
*   query<br/>
    A parallelized BLAST wrapper that allows you to run a set of query sequences (nucleotide or protein) against one or many databases in batch.
    Includes:
    *   Support for blastn, tblastn, and blastp based on input types (auto-detected).
    *   Optional reporting of only strongest match per query per genome.
    *   Automatic filtering based on identity, coverage and/or e-value thresholds.
    *   Output includes both full and filtered results tables, plus alignment files if desired.
    *   Motif mining for amino acid queries.
    *   Visualization of gene hits and sequence matches (hi-res PNG or PDF)
*   Motif support in query: <br/>
    When using blastp (amino acid query against protein database), users may specify one or more amino acid motif(s) (e.g., WXWXIP (single motif) | XAXH GHXXGE (multiple motifs, space-separated list)) using the `--motif` flag. This performs a regex-based search across all BLAST hits, independent from internally-curated or user-defined percent identity, query coverage, and e-value thresholds, ensuring detection of conserved motifs even in low-identity or heterologous alignments that might otherwise be filtered out. This is particularly useful for detecting signature domains (e.g., catalytic triads, DNA-binding motifs) in diverse sequence families. Whole query matches or just the motif string may be exported to FASTA using `--motif-fasta-out` with or without `--motif-only`.
    *    For focused motif investigations across multiple queries, each motif's associated query file may be linked in the command string. This will result in only results associated with that query file being parsed for motif matches.
    *    Example:
            *    Query files: AT_domain.faa, KS_domain.faa
            *    Motifs: GHXXGE, TAXXSS
            *    Unlinked Query motif search: <br/>
                `seqforge query -d <DB_DIR> -q <Query_DIR> -o <RESULTS_DIR> -f <FASTA_DIR> --motif GHXXGE TAXXSS`
                    * This results in all queries being searched for each motif
            *    Linked Query motif search: <br/>
                `seqforge query -d <DB_DIR> -q <Query_DIR> -o <RESULTS_DIR> -f <FASTA_DIR> --motif GHXXGE{AT_domain} TAXXSS{KS_domain}` <br/>
                    *    This results in only AT_domain.faa being searched for GHXXGE and only KS_domain.faa being searched for TAXXSS.
                    *    The exact query file base name should be placed within brackets as a part of the motif string. 
*   BLAST hits and motif matches may be visualized using the `--visualize` flag. BLAST hits across genomes are represented by a heatmap, where the color intensity of each individual cell reflects that query's percent identity within a specific genome. For queries that return > 1 hit per genome, only  the strongest hit will be used to construct the heatmap. If `--visualize` is combined with `--motif`, motif matches will be illustrated as a sequence logo representing amino acid frequencies. Plots may be saved as either a high resolution PNG or in PDF format (recommended). 
______________________________________________________________________________________________________________________________________
Example Plots:

<img width="21240" height="4457" alt="query_identity_heatmap" src="https://github.com/user-attachments/assets/d07c9521-c99f-4ddf-8f14-6e43f92d7340" />

<img width="2970" height="716" alt="motif_1_logo" src="https://github.com/user-attachments/assets/92142974-6c20-4bb1-8c20-b7be5951787f" />

______________________________________________________________________________________________________________________________________
### <ins>Module 2: Sequence Investigation<ins/>

Purpose: Extract meaningful biological context from BLAST hits <br/>
*   extract <br/>
    Extract aligned sequences identified via the Genome Search Query pipeline from original FASTA files.
    Features:
    *   Optional translation of nucleotide hits to protein for full gene alignments.
        *    This should only be used with full gene alignments.
    *   Upstream/downstream padding for context-based analysis.
    *   Filtering using percent identity, query coverage, and/or e-value.
*   extract-contig<br/>
    Extract *entire contigs* from reference assemblies based on where BLAST hits occurred. Ideal for identifying genomic context of hits in metagenomic assemblies too large to open via a genome browser.
    *   For .faa files containing coding sequences, `extract-contig` will extract the full protein sequence from the source file, not just the aligned region identified via `query` for query coverage < 100%. 
______________________________________________________________________________________________________________________________________
### <ins>Module 3: Utilities<ins/>

Purpose: General-use tools for various genomic workflows <br/>
*   sanitize<br/>
    *   Remove special characters from input files—**required** for all Genome Search analyses.
    *   Choose to rename files using `--in-place` (strongly recommended).
    *   Or copy files to a new directory and rename the copies, leaving the original files untouched.
*   split-fasta<br/>
    Split multi-FASTA files into smaller chunks for downstream processing.
    *   Choose fixed sequence count per chunk or split one sequence per file.
    *   Optional compress output for storage/transfer efficiency.
*   search<br/>
    Extract isolation metadata from **GenBank** or **JSON** files.
    *   User-defined or comprehensive field extraction to CSV or TSV.
*   fasta-metrics<br/>
    Compute common assembly metrics from an input FASTA file or all FASTA files within a directory
    *   The following metrics are calculated
        *   Genome size (bp)
        *   Number of contigs
        *   Longest contig
        *   Shortest contig
        *   Basic contig size distribution
        *   Total length of contigs ≥ 1 kb
        *   GC content (%)
        *   N count
        *   auN
        *   N50, N90
        *   L50, L90
        *   Lengths of each contig (CSV only)
* unique-headers <br/>
    Append unique identifiers to each FASTA header within a multi-FASTA file
    *   Submit a single FASTA file or a directory of FASTA files
    *   Appends source file name + unique alphanumeric barcode to header
    *   Example:
        *   Source file: GCA_000346795_1.faa
        *   Example header: >hypothetical
        *   unique-headers modification: >hypothetical_GCA_000346795_1_54uMe
______________________________________________________________________________________________________________________________________
______________________________________________________________________________________________________________________________________

## Recommended Workflow:

```mermaid
graph TD;
    Z{{FASTA}}-->A{Sanitize};
    A{Sanitize}-->B{{Database-Creation}};
    B{{Database-Creation}}-->C[Query];
    C[Query]-->D[Extract];
    C[Query]-->E[Extract-Contig];
    A{Sanitize}-->F[FASTA-Metrics];
    A{Sanitize}-->H[Split-FASTA];
    I[JSON/GenBank]-->J[Search];

    %% Custom style for Sanitize node
    style Z fill:#11d393,stroke:#fff5ee,stroke-width:4,font-size:30,color:#454545
    style A fill:#d48074,stroke:#fff5ee,stroke-width:3,font-size:26px
    style B fill:#d40078,stroke:#FCF5E5,stroke-width:2,font-size:22px
    style C fill:#920075,stroke:#333,stroke-width:1,font-size:18px
    style D fill:#650D89,stroke:#333,stroke-width:1,font-size:18px
    style E fill:#023788,stroke:#333,stroke-width:1,font-size:18px
    style F fill:#fd1d53,stroke:#333,stroke-width:1,font-size:18px
    style H fill:#2e2157,stroke:#333,stroke-width:1,font-size:18px
    style I fill:#11c9d3,stroke:#333,stroke-width:1,font-size:22px,color:#454545
    style J fill:#6e1515,stroke:#333,stroke-width:1,font-size:18px
```
______________________________________________________________________________________________________________________________________
**All pipelines** should start with a check for special characters in input FASTA file names. Special characters are defined as any 
character aside from a-z, A-Z, 0-9 and underscores. The BLAST+ architecture cannot parse special characters; inclusion of these
characters in a filename (e.g., GCA_900638215.1.fna, the **non-extension period** is pipeline-breaking) will result in database creation failure.

To scan your FASTA population, simply run: `seqforge --sanitize -f /path/to/FASTA/file(s) -e fasta --dry-run`.
This will print any problematic filenames to the console without making changes. If changes are needed, re-run without `--dry-run` 
using either `--in-place` (recommended) or `--sanitize-outdir <dir>`. See Module 3: Utilities > seqforge sanitize for more details.

Periods, hyphens, colons, semicolons, and whitespace are replaced by '_'. Parenthesis and quotation marks are deleted. GCA_900638215.1.fna
(**non-compliant**) becomes GCA_900638215_1.fna (**compliant**).
______________________________________________________________________________________________________________________________________
**Suggested Workspace:** <br/>
For organization, we suggest generating individual directories for each data type prior to the start of any SeqForge module:

`cd /path/to/working/directory` <br/>
`mkdir -p FASTA DBs Query Results`, **where:**

**DBs:** Output directory to store BLAST databases generated via `seqforge makedb` <br/>
**FASTA:** Directory of genomes to parse, used to create BLAST databases <br/>
**Query:** Directory containing query file(s) to use with `seqforge query` <br/>
**Results:** Directory where Genome Search results will be stored

We additionally recommend generating output directories for any additional SeqForge modules used. 
______________________________________________________________________________________________________________________________________
______________________________________________________________________________________________________________________________________

**Author:** <br />
    Elijah R. Bring Horvath, PhD (https://github.com/ERBringHorvath)

**License:** <br />
    This program is shared under MIT License, which allows for modification and redistribution with attribution.
______________________________________________________________________________________________________________________________________

# Installation

### Install NCBI BLAST+

**SeqForge uses [NCBI BLAST+](https://pubmed.ncbi.nlm.nih.gov/20003500/)**

Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, **2009**. <br />
BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. doi:10.1186/1471-2105-10-421

**Install BLAST+**

Download latest version of [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

Or using Conda:

1. Install Conda [miniforge](https://github.com/conda-forge/miniforge/) if not already installed

2. Create Conda environment

`conda create -n seqforge`

3. Activate Conda envrionment

`source activate seqforge`

4. Install BLAST+

`conda -y install bioconda::blast`

**Verify BLAST Installation**

`makeblastdb -h` <br />
`blastn -h`

If these commands run without error, BLAST is correctly installed. If an error occurs, refer to the [BLAST+ documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/index.html#index)

## SeqForge Installation

If not already installed, install [Git](https://github.com/git-guides/install-git) <br />
* Linux/macOS systems should have this installed by default <br />
* To test installation, open the terminal and type `git --version` <br />
* For macOS users, you should see something like `git version 2.37.1 (Apple Git-137.1)`

### **Clone SeqForge from source**

We suggest installing SeqForge within your Home folder, such as `/home/user/` 

Change directory to desired installation path

`cd /home/user`

Clone SeqForge from the repository

`git clone https://github.com/ERBringHorvath/SeqForge`

Add SeqForge to your PATH

1. Open your profile in a text editor. This might be `~/.bash_profile`, `~/.bashrc`, or `~/.zshrc`
2. Add the following line to the end of the file:

`export PATH=$PATH:/home/user/SeqForge/seqforge`

Replace `/home/user/SeqForge/seqforge` with the actual path to the directory containing the executable. <br />
Whatever the initial directory, this path should end with `/SeqForge/seqforge`

Save the file and restart your terminal or run `source ~/.bashrc` or `~/.bash_profile` (Linux/macOS) or `source ~/.zshrc` (macOS)

**Install Dependencies**

`pip install -r requirements.txt` or `pip3 install -r requirements.txt`

**Verify SeqForge Installation**

`seqforge --module-health` <br/>

You should see: <br/>
```
Module Status Report:
 - makedb: Available
 - query: Available
 - split_fasta: Available
 - extract_sequences: Available
 - extract_contigs: Available
 - mask: Available
 - search: Available
 - sanitize: Available
 - fasta_metrics: Available
 - generate_unqiue_fasta_headers: Available
```
If there is a problem, the status will read 'Broken or Missing' <br/>
Submit issues [here](https://github.com/ERBringHorvath/SeqForge/issues)

Help menu and version: <br/>
`seqforge --help` <br />
`seqforge --version`

Indifivual module help: <br/>
`seqforge <module> -h`

NOTE: Permissions should automatically be applied during installation. If you get a `permission denied` message when running `seqforge`
permissions may need to be changed manually. To do this, you can use the following command (requires administrator privileges):

`chmod +x /path/to/SeqForge/seqforge/seqforge`

______________________________________________________________________________________________________________________________________
______________________________________________________________________________________________________________________________________

# <ins>Module 1: Genome Search</ins>

**Building a BLAST+ Database Library**

seqforge makedb: <br/>
**Required arguments:** <br/>
`-f`, `--fasta-directory`: path to the directory containing input files in FASTA format <br />
`-o`, `--out`: path to directory where you want to store your databases

**Optional arguments:** <br/>
`-T`, `--threads`: number of cores to dedicate for multiprocessing (default = 4) <br/>
`-s`, `--sanitize`: remove pipeline-breaking special characters from file names (renames in-place) <br/>
`--temp-dir`: specify a temporary directory for archive extraction (default = /tmp/) <br/>
`--keep-temp-files`: for archive submission; retains temporary directory generated at /tmp/seqforge_fasta_extract_* <br/>
`-p`, `--progress`: progress reporting mode; choices = 'bar', 'verbose', 'none'. --progress verbose prints a line per item.

Example: <br />
`seqforge makedb -f /path/to/FASTA/files -o /path/to/results/folder -T 8`

**Database type is automatically detected during database creation following standard FASTA extension practices** <br/>
.fasta, .fa, .fas, .ffn, .fna == nucleotide <br/>
.faa == protein 
______________________________________________________________________________________________________________________________________

**Querying a database library**

seqforge query: <br />
**Required arguments:** <br/>
`-d`, `--database`: path to directory containing BLAST+ databases <br />
`-q`, `--query-files`: path to directory containing query files in amino acid FASTA format <br />
`-o`, `--output`: path to directory to store results

**Required for nucleotide queries (blastn):** <br/>
`-N`, `--nucleotide-query`: use blastn for queries in nucleotide FASTA format

**Optional arguments:** <br/>
`-T`, `--threads`: number of cores to dedicate for multiprocessing (default = 4) <br />
`-R`, `--report-strongest-match`: report only the single strongest match for each query <br />
`--min-perc`: define minimum percent identity threshold. Default = 90 <br />
`--min-cov`: define minimum query coverage threshold. Default = 75 <br />
`-e`, `--evalue`: maximum e-value cutoff, default 0.00001 <br />
`--min-seq-len`: define minimum sequence length for short nucleotide sequence queries (use with caution and only with `--nucleotide-query`) <br />
`-a`, `--no-alignment-files`: suppress alignment file creation <br/>
`--temp-dir`: specify a temporary directory for archive extraction (default = /tmp/) <br/>
`--keep-temp-files`: retain individual *_results.txt files in output directory <br/>
`--motif`: amino acid motif (e.g., WXWXIP or space-separated list) to search within blastp hits. X is treated as a wildcard. Only for use with blastp queries <br/>
`--motif-fasta-out`: export motif query alignments to FASTA <br/>
`--motif-only`: for use with `--motif-fasta-out`; export only motif string to FASTA <br/>
`-f`, `--fasta-directory`: path to FASTA file(s) used to create BLAST databases. Required if using `--motif` <br/>
`-V`, `--visualize`: generate heatmap of BLAST hits and sequence logo of motif hits if `--motif` returns matches <br/>
`-P`, `--pdf`: override PNG output of visualize and instead generated a PDF (use in combination with `--visualize`) <br/>
`-p`, `--progress`: progress reporting mode; choices = 'bar', 'verbose', 'none'. --progress verbose prints a line per item.

Basic example: <br />
`seqforge query -T 8 -d /path/to/blast/database/files -q /path/to/query/files/` <br /> 
`-o /path/to/results/folder`

**Witches' Brew:** <br/>
`seqforge query -T 8 \` <br/>
`-d /path/to/blast/databases \` <br/>
`-q /path/to/query/files \` <br/>
`-o /path/to/results/directory \` <br/>
`--min-perc 75` `--min-cov 70 --evalue 0.001 \` <br/>
`--report-strongest-matches \` <br/>
`--motif TAXXSS{query1} GHXXGE{query2} -f /path/to/FASTA/files \` <br/>
`--motif-fasta-out --visualize --pdf \` <br/>
`--motif-fasta-out`

All Query results are concatenated to `all_results.csv` and either `all_filtered_results.csv` or `filtered_results.csv`. 
All plots and files are saved to the output directory specified by `--output`.
______________________________________________________________________________________________________________________________________
______________________________________________________________________________________________________________________________________

# <ins>Module 2: Sequence Investigation</ins>

## Extract Sequences from a SeqForge Query

seqforge extract: <br />
**Required arguments:** <br/>
`-c`, `--csv-path`: path to results csv file from `seqforge query` <br />
`-f`, `--fasta-directory`: path to reference FASTA assemblies <br />
These should be the FASTA files the BLAST databases were created from and should have the same basename as the query results files <br />
`-o`, `--output-fasta`: output file to contain sequences, defaults to current working directory

**Optional arguments:** <br/>
`-T`, `--threads`: number of cores to dedicate for multiprocessing (default = 1) <br />
`-e`, `--evalue`: maximum e-value threshold, default = 0.00001 <br />
`--min-perc`: minimum percent identity threshold. Default = 90 <br />
`--min-cov`: minimum query coverage threshold. Default = 75 <br />
`--translate`: translates extracted nucleotide sequence(s) <br/>
`--up`: extract additional basepairs upstream of aligned sequence <br/>
`--down`: extract additional basepairs downstream of aligned sequence <br/>
`--temp-dir`: specify a temporary directory for archive extraction (default = /tmp/) <br/>
`--keep-temp-files`: for archive submission; retains temporary directory generated at /tmp/seqforge_fasta_extract_*

**NOTE:** Translation of sequences is optional, however care should be used when translating extracted nucleotide sequences, as BLAST results may not always contain a full CDS. To allow for this, when the `--translate` argument is called, extracted sequences will be trimmed to only include complete codons, which may affect interpretation of results. Generally, BLAST+ will return sequence coordinates in the correct frame; this has been tested and no mistranslations have been logged, but there are always exceptions. 

**NOTE:** `--up` and `--down` flags are incompatible with `--translate`, as 6-frame translation is not currenlty supported.

**NOTE:** Results files and FASTA reference assemblies <ins>**must**</ins> share the same basename:

Example basename: 'FILE' <br />
Example FASTA: FILE.fasta <br />
Example results file: FILE_results.txt

If SeqForge is used for database creation and queries, matching basenames will be generated automatically

**Example usage:** <br />
`seqforge extract -c /path/to/results/file -f /path/to/reference/FASTA/files -T 8 -o sequences.fa`

**Witches' Brew:** <br/>
`seqforge extract -T 8 \` <br/>
`-c /path/to/results/file -f /path/to/reference/FASTA/files \` <br/>
`-o sequences.fa` <br/>
`--min-perc 75 --min-cov 70 --evalue 0.001 \` <br/>
`--up 1200 --down 1200`

**NOTE:** <br />
Results file needs to be `all_results.csv`, `all_filtered_results.csv`, or `filtered_results.csv`, which are automatically generated using `seqforge query`

If percent identity and query coverage were set manually during `seqforge query`, these values will need to be reflected when using `mutliblast extract` using `--min-perc` and/or `--min-cov` <br />
For instance, if `seqforge query` was called using `--perc 75`, but the `seqforge extract` minimum percent identity is left at its default value (90), the appropriate sequences may not be extracted, as they may fall beneath the internally curated `--min-perc` theshold. 

`seqforge extract` will generate a multi-FASTA file of all sequences identified by `seqforge queryP`/`query` based on the default or user-defined e-value cutoff.
______________________________________________________________________________________________________________________________________

## Extract Entire Contig ##

seqforge extract-contig: <br />
**Required arguments:** <br/>
`-c`, `--csv-path`: path to csv results file from `seqforge query` <br />
`-f`, `--fasta-directory`: path to reference FASTA assemblies <br />
These should be the FASTA files the BLAST databases were created from and should have the same basename as the query results files <br />
`-o`, `--output-fasta`: output file to contain sequences, defaults to current working directory

**Optional arguments:** <br/>
`-T`, `--threads`: number of cores to dedicate, default is 1 <br />
`-e`, `--evalue`: maximum e-value threshold, default = 0.00001 <br />
`--min-perc`: minimum percent identity threshold. Default = 90 <br />
`--min-cov`: minimum query coverage threshold. Default = 75 <br />
`--temp-dir`: specify a temporary directory for archive extraction (default = /tmp/) <br/>
`--keep-temp-files`: for archive submission; retains temporary directory generated at /tmp/seqforge_fasta_extract_*

**NOTE:** Results files and FASTA reference assemblies <ins>**must**</ins> share the same basename for both `seqforge extract` and `seqforge extract-contig`:

Example basename: 'FILE' <br />
Example FASTA: FILE.fasta <br />
Example results file: FILE_results.txt

If SeqForge is used for database creation and queries, matching basenames are handled automatically

**Example usage:** <br />
`seqforge extract-contig -d /path/to/results/files -f /path/to/reference/FASTA/files -T 8 -o contigs.fa`

`seqforge extract-contig` will generate a multi-FASTA file of all contigs harboring a matching <br /> 
sequence identified by `seqforge query` based on the default or user-defined thresholds. 

This program was designed for use with metagenome mining, as metagenomic assemblies are often too large to explore using a genome browser. If short-read assembly methods are used, contigs harboring genes of interest may be extracted; contigs will likely be more tractible to parsing using a genome browser if manual annotation is needed. 

______________________________________________________________________________________________________________________________________
______________________________________________________________________________________________________________________________________

# <ins>Module 3: Utilties<ins/>

## Sanitize File Names

seqforge sanitize: <br/>
`-i`, `--input`: FASTA file or directory of FASTA files to be cleaned <br/>
`-e`, `--extension`: extention of files to sanitize; for all standard FASTA formats, simply pass `-e fasta` <br/>
`-I`, `--in-place`: overwrite problematic file names with names containing no special  characters (recommended) <br/>
`-S`, `--sanitize-outdir`: path to new directory for files. This option copies file(s) to a new directory and removes special characters from
        the copied file names, leaving the original file names intact (storage intensive) <br/>
`--dry-run`: preview changes without committing <br/>
______________________________________________________________________________________________________________________________________

## FASTA File Metrics

seqforge fasta-metrics: <br/>
`-f`, `--fasta-directory`: path to FASTA file or directory of FASTA files to be analyzed <br/>
`-o`, `--output`: optional name for CSV summary (default: fasta_metrics_summary.csv) <br/>
`-M`, `--min-contig-size`: minimum contig size to include for all calculations (default = 500) <br/>
`--temp-dir`: specify a temporary directory for archive extraction (default = /tmp/) <br/>
`--keep-temp-files`: for archive submission; retains temporary directory generated at /tmp/seqforge_fasta_extract_*
______________________________________________________________________________________________________________________________________

## Split Multi-FASTA files

seqforge split-fasta: <br />
`-i`, `--input`: input multi-FASTA file <br />
`-o`, `--output-dir`: output directory for split FASTA files <br />
`-F`, `--fragment`: split multi-FASTA file into defined chunks of *n* sequences each <br />
`-C`, `--compress`: compress output files as .gz <br />
______________________________________________________________________________________________________________________________________

## Extract Sequence Metadata from JSON/GenBank Files

seqforge search: <br/>
`-i`, `--input`: input file(s) (.json or .gb/.gbk)<br/>
`-o`, `--output`: output file (e.g., .csv, .tsv, .json) <br/>
`--json`: parse only JSON files in an input directory <br/>
`--gb`: parse only GenBank files in an input directory <br/>
`--all`: extract all available metadata <br/>
`--fields`: space-separated list of metadata fields to extract

**Available fields:** <br/>
    accession, organism, strain, isolation_source, host, region, lat_lon, collection_date, <br/> 
    collected_by, tax_id, comment, keywords, sequencing_tech, release_date

**Example usage:** <br/>
`seqforge search -i /path/to/input.(json | gbk) -o metadata.csv --fields accession isolation_source host region` <br/>
`seqforge search -i /path/to/input/files -o metadata.csv --all --json`
______________________________________________________________________________________________________________________________________

## Generate Unique FASTA Headers 
### (AUGUSTUS base prediction compliant)

seqforge unique-headers: <br/>
`-i`, `--input`: path to FASTA file(s) <br/>
`-o`, `--output-dir`: directory for output FASTA files (unless using `--in-place`) <br/>
`-I`, `--in-place`: modify input files in-place (uses temporary files for safety) <br/>
`-p`, `--progress`: progress reporting mode; choices = 'bar', 'verbose', 'none'. --progress verbose prints a line per record.
______________________________________________________________________________________________________________________________________
______________________________________________________________________________________________________________________________________

# Citations

Cite SeqForge: <br />
SeqForge (https://github.com/ERBringHorvath/SeqForge)

Cite NCBI BLAST+: <br />
Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, **2009**. <br />
BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. doi:10.1186/1471-2105-10-421

Cite Logomaker: <br/>
Ammar Tareen, Justin B Kinney, **2020** <br/>
Logomaker: beautiful sequence logos in Python. *Bioinformatics*, 36, 7, 2272–2274, doi.org/10.1093/bioinformatics/btz921
