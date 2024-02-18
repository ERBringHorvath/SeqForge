# **multiBLAST**

# Install NCBI BLAST+

*multiBLAST uses NCBI BLAST+ Software*

Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, 2009. 
BLAST+: architecture and applications. BMC Bioinformatics, 10, 421. doi:10.1186/1471-2105-10-421

**Install BLAST+**

https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Or using Conda

1. Create Conda environment

`conda create -n blast`

2. Activate Conda envrionment

`source activate blast`

3. Install BLAST+

`conda -y install bioconda::blast`

**Verify BLAST Installation**

`makeblastdb -h` \n
`blastn` -h

# multiBLAST Installation

We suggest installing multiBLAST within your Home folder, such as `/Users/your_username/multiblast` 

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

**Verify 
