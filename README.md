# **multiBLAST**

**This program automates the process of running BLAST (Basic Local Alignment Search Tool) <br />
queries against multiple databases and organizing the results. It is designed to handle various <br />
types of BLAST searches (e.g., blastn, tblastn) and formats the output for easy analysis.** <br />

# Install NCBI BLAST+

*multiBLAST uses [NCBI BLAST+ Software](https://pubmed.ncbi.nlm.nih.gov/20003500/)*

Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, 2009. <br />
BLAST+: architecture and applications. BMC Bioinformatics, 10, 421. doi:10.1186/1471-2105-10-421

**Install BLAST+**

Download latest version of [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

Or using Conda

1. Create Conda environment:

`conda create -n blast`

2. Activate Conda envrionment:

`source activate blast`

3. Install BLAST+:

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

`multiblast --help`


