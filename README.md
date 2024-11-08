# **<ins>multiBLAST<ins/>**

**multiBLAST automates the process of running multiple BLAST (Basic Local Alignment Search Tool) queries against multiple
databases and organizing the results. It is designed to handle various types of BLAST searches (nucelotide, translated nucleotide)
and formats the output for easy analysis. Additionally, aligned sequences may be extracted from genome assemblies. To facilitate
metagenomic analyses, contigs harboring genes of interest identified via multiBLAST may additionally be extracted for further analysis.
multiBLAST is scalable and has been used to run dozens of gene queries against hundreds of genomes.** 

**Author:** <br />
    Elijah R. Bring Horvath, PhD (https://github.com/ERBringHorvath)

**License:** <br />
    This script is shared under MIT License, which allows for modification and redistribution with attribution.

## Install NCBI BLAST+

**multiBLAST uses [NCBI BLAST+](https://pubmed.ncbi.nlm.nih.gov/20003500/)**

Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, **2009**. <br />
BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. doi:10.1186/1471-2105-10-421

**Install BLAST+**

Download latest version of [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

Or using Conda:

1. Install Conda [miniforge](https://github.com/conda-forge/miniforge/) if not already installed

2. Create Conda environment

`conda create -n blast`

3. Activate Conda envrionment

`source activate blast`

4. Install BLAST+

`conda -y install bioconda::blast`

**Verify BLAST Installation**

`makeblastdb -h` <br />
`blastn -h`

If these commands run without error, BLAST is correctly installed. If an error occurs, refer to the [BLAST+ documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/index.html#index)

## multiBLAST Installation

If not already installed, install [Git](https://github.com/git-guides/install-git) <br />
* Linux/Unix systems should have this installed by default <br />
* To test installation, open the terminal and type `git --version` <br />
* For macOS users, you should see something like `git version 2.37.1 (Apple Git-137.1)`

## **Clone multiBLAST from source**

We suggest installing multiBLAST within your Home folder, such as `/Users/user/` 

Change directory to desired installation path

`cd /Users/user`

Clone multiBLAST from the repository

`git clone https://github.com/ERBringHorvath/multiBLAST`

Add multiBLAST to your PATH

1. Open your profile in a text editor. This might be `~/.bash_profile` or `~/.zshrc`
2. Add the following line to the end of the file:

`export PATH=$PATH:/Users/user/multiblast/bin`

Replace `/Users/user/multiblast/bin` with the actual path to the directory containing the executable. <br />
Whatever the initial directory, this path should end with `/multiblast/bin`

Save the file and restart your terminal or run `source ~/.bash_profile` (Linux/Unix) or `source ~/.zshrc` (macOS)

**Install Dependencies**

`pip install -r requirements.txt` or `pip3 install -r requirements.txt`

**Verify multiBLAST Installation**

`multiblast --help` <br />
`multiblast --version`

NOTE: Permissions should automatically be applied during installation. If you get a `permission denied` message when running `multiblast` <br />
permissions may need to be changed manually. To do this, you can use the following command:

`chmod +x /path/to/multiblast/bin/multiblast`

## Example Usage

**Building a BLAST+ Database Library**

multiblast makedb: <br />
`-f`, `--file_directory`: path to the directory containing input files in FASTA format <br />
`-d`, `--dbtype`: specify what sort of database you want to create (`nucl`, nucleotide, `prot`, protein) <br />
`-o`, `--out`: path to directory where you want to store your databases

Example: <br />
`multiblast makedb -f /path/to/FASTA/files -d nucl -o /path/to/results/folder`

**Querying a database library**

multiblast query: <br />
`-d`, `--database`: path to directory containing BLAST+ databases <br />
`-q`, `--query_files`: path to directory containing query files in amino acid FASTA format <br />
`-e`, `--evalue`: maximum e-value cutoff, default 0.00001 <br />
`-o`, `--output`: path to directory to store results
`-T`, `--threads`: number of cores to dedicate for multi-threading <br />
`--report-strongest-match`: report only the single strongest match for each query <br />
`--perc`: define minimum percent identity threshold. Default = 90 <br />
`--cov`: define minimum query coverage threshold. Default = 75 <br />
`--nucleotide-query`: use blastn for queries in nucleotide FASTA format <br />
`--min-seq-len`: define minimum sequence length for short amino acid/nucleotide sequence queries (use with caution) <br />

Example: <br />
`multiblast query -T 8 -d /path/to/blast/database/files -q /path/to/query/files/` <br /> 
`-o /path/to/results/folder`

All multiBLAST results are concatenated to `all_results.csv` and either `all_filtered_results.csv` or <br /> 
`filtered_results.csv` within the output folder designated by `-o, --output`

# <ins>Accessory Scripts</ins>

## Extract Sequences from a multiBLAST Query

multiblast extract: <br />
`-d`, `--results_directory`: path to directory containing multiBLAST results files <br />
`-f`, `--fasta_directory`: path to reference FASTA assemblies <br />
&nbsp;&nbsp;&nbsp;&nbsp;These should be the FASTA files the BLAST databases were created from and should have the same basename as the query results files <br />
`-o`, `--output_fasta`: output file to contain sequences, defaults to current working directory <br />
`-T`, `--threads`: number of cores to dedicate, default is 1 <br />
`--min-evalue`: maximum e-value threshold, default = 0.00001 <br />
`--min-perc`: minimum percent identity threshold. Default = 90 <br />
`--min-cov`: minimum query coverage threshold. Default = 75 <br />
`--translate`: translates extracted nucleotide sequence(s)

**NOTE:** Translation of sequences is optional, however care should be used when translating extracted nucleotide sequences, as BLAST results may not always contain a full CDS. To allow for this, when the `--translate` argument is called, extracted sequences will be trimmed to only include complete codons, which may affect interpretation of results.

**NOTE:** Results files and FASTA reference assemblies <ins>**must**</ins> share the same basename:

Example basename: 'FILE' <br />
&nbsp;&nbsp;&nbsp;&nbsp;Example FASTA: FILE.fasta <br />
&nbsp;&nbsp;&nbsp;&nbsp;Example results file: FILE_results.txt

If multiBLAST is used for database creation and queries, matching basenames should be generated automatically

**Example usage:** <br />
`multiblast extract -d /path/to/results/file -f /path/to/reference/FASTA/files -T 8 -o sequences.fa`

**NOTE:** <br />
Results file should be `all_results.csv`, `all_filtered_results.csv`, or `filtered_results.csv`, which are automatically generated using `multiblast query`

If percent identity and query coverage were set manually during `multiblast query`, these values will need to be reflected when using `mutliblast extract` using `--min-perc` and/or `--min-cov` <br />
For instance, if `multiblast query` was called using `--perc 75`, but the `multiblast extract` minimum percent identity is left at its default value (90), the appropriate sequences may not be extracted, as they may fall beneath the internally curated `--min-perc` theshold. 

`multiblast extract` will generate a multi-FASTA file of all sequences identified by `multiblast queryP`/`query` based on the default or user-defined e-value cutoff.

## Extract Entire Contig ##

`multiblast extract-contig`: <br />
`-d`, `--results_directory`: path to directory containing multiBLAST results files <br />
`-f`, `--fasta_directory`: path to reference FASTA assemblies <br />
&nbsp;&nbsp;&nbsp;&nbsp;These should be the FASTA files the BLAST databases were created from and should have the same basename as the query results files <br />
`-o`, `--output_fasta`: output file to contain sequences, defaults to current working directory <br />
`-T`, `--threads`: number of cores to dedicate, default is 1 <br />
`--min-evalue`: maximum e-value threshold, default = 0.00001 <br />
`--min-perc`: minimum percent identity threshold. Default = 90 <br />
`--min-cov`: minimum query coverage threshold. Default = 75 <br />

**NOTE:** Results files and FASTA reference assemblies <ins>**must**</ins> share the same basename for both `multiblast extract` and `multiblast extract-contig`:

Example basename: 'FILE' <br />
&nbsp;&nbsp;&nbsp;&nbsp;Example FASTA: FILE.fasta <br />
&nbsp;&nbsp;&nbsp;&nbsp;Example results file: FILE_results.txt

If multiBLAST is used for database creation and queries, matching basenames are handled automatically

**Example usage:** <br />
`multiblast extract-contig -d /path/to/results/files -f /path/to/reference/FASTA/files -T 8 -o contigs.fa`

`multiblast extract-contig` will generate a multi-FASTA file of all contigs harboring a matching <br /> 
sequence identified by `multiblast query` based on the default or user-defined thresholds. 

This program was designed for use with metagenome mining, as metagenomic assemblies are often too large to explore using a genome browser. If short-read assembly methods are used, contigs harboring genes of interest may be extracted; contigs will likely be more tractible to parsing using a genome browser if manual annotation is needed. 

## Split multi-FASTA files

Often times when downloading large genomic datasets, individual FASTA files will be concatenated into one large multi-FASTA file. This script is designed to split multi-FASTA files into individual FASTA files for use with the multiBLAST or other bioinformatic platforms.

**Example usage:**

`multiblast split_fasta -i /path/to/multiFASTA/file -o /path/to/results/folder`

# Citations

Cite multiBLAST: <br />
multiBLAST (https://github.com/ERBringHorvath/multiBLAST)

Cite NCBI BLAST+: <br />
Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, **2009**. <br />
BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. doi:10.1186/1471-2105-10-421

