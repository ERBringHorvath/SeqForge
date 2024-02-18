# **multiBLAST**

**This program automates the process of running BLAST (Basic Local Alignment Search Tool) queries against multiple <br />
databases and organizing the results. It is designed to handle various types of BLAST searches (e.g., blastn, tblastn) <br />
and formats the output for easy analysis.** 

**Author:** <br />
    Elijah R. Bring Horvath (https://github.com/ERBringHorvath)

**License:** <br />
    This script is shared under MIT License, which allows for modification and redistribution with attribution.

**Note**: <br />
    This script is intended for research and academic purposes. 
    Please ensure you have the necessary permissions to use the databases and query files with BLAST.

# Install NCBI BLAST+

*multiBLAST uses [NCBI BLAST+ Software](https://pubmed.ncbi.nlm.nih.gov/20003500/)*

Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, 2009. <br />
BLAST+: architecture and applications. BMC Bioinformatics, 10, 421. doi:10.1186/1471-2105-10-421

**Install BLAST+**

Download latest version of [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

Or using Conda:

1. Create Conda environment

`conda create -n blast`

2. Activate Conda envrionment

`source activate blast`

3. Install BLAST+

`conda -y install bioconda::blast`

**Verify BLAST Installation**

`makeblastdb -h` <br />
`blastn` -h

If these commands run without error, BLAST is correctly installed. If an error occurs, refer to the [BLAST+ documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/index.html#index)

# multiBLAST Installation

We suggest installing multiBLAST within your Home folder, such as `/Users/your_username/` 

Change directory to desired installation path

`cd /Users/your_username`

Clone multiBLAST from the repository

`git clone https://github.com/ERBringHorvath/multiBLAST`

Add multiBLAST to your PATH

1. Open your profile in a text editor. This might be `~/.bash_profile` or `~/.zshrc`
2. Add the following line to the end of the file:

`export PATH=$PATH:/path/to/multiblast`

Replace `/path/to/multiblast` with the actual path to the directory containing executable

Save the file and restart your terminal or run `source ~/.bash_profile` (Linux/Unix) or `source ~/.zshrc` (macOS)

**Install Dependencies**

`pip install -r requirements.txt` or `pip3 install -r requirements.txt`

**Verify multiBLAST Installation**

`multiblast --help` <br />
`multiblast --version`

# Example Usage

**Building a BlAST+ Database Library**

multiblast makedb: <br />
`-f, --file_directory`: path to the directory containing input files in FASTA format <br />
`-d, --dbtype`: specify what sort of database you want to create (`nucl`, nucleotide, `prot`, protein) <br />
`-o, --out`: path to directory where you want to store your databases

Example: <br />
`multiblast makedb -f /path/to/FASTA/files/folder -d nucl -o /path/to/results/folder`

**Querying a database library**

multiblast query: <br />
`-m, --method`: BLAST method to perform <br />
    `tblastn`, search protein query through nucleotide database <br />
    `blastn`, search nucleotide query trhough nucleotide database <br />
`-d, --database`: path to directory containing BLAST+ databases <br />
`-q, --query_files`: path to directory containing query files in FASTA format <br />
`-e, --evalue`: maximum e-value cutoff <br />
`-o, --output`: path to directory to store results

Example: <br />
`multiblast query -m tblastn -d /path/to/blast/database/folder -q /path/to/query/files/folder -e 0.01 -o /path/to/results/folder`

All multiBLAST results are concatenated to `multiblast_results.csv` within the output folder designated by `-o, --output`

# Citations

Cite multiBLAST: <br />
multiBLAST v.1.0.0 (https://github.com/ERBringHorvath/multiBLAST)

Cite NCBI BLAST+: <br />
Camancho C et al., **2009**. BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. 

