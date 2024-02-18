# multiBLAST
"""
multiBLAST

This program automates the process of running BLAST (Basic Local Alignment Search Tool) 
    queries against multiple databases and organizing the results. It is designed to handle 
    various types of BLAST searches (e.g., blastn, tblastn) and formats the output for easy analysis.

multiBLAST utilizes NCBI BLAST+: 
    Camancho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL, 2009. 
    BLAST+: architecture and applications. BMC Bioinformatics, 10, 421. doi:10.1186/1471-2105-10-421

Usage:
    - Ensure BLAST+ command line tools are installed and accessible in your environment.
    - Prepare your query files and BLAST database directories.
    - Run the script and follow the prompts to specify the BLAST method, paths to query files, 
        database directories, and the desired output location.

Requirements:
    - Python 3.x
    - pandas library
    - BLAST+ command line tools installed and in your system's PATH.

Features:
    - Supports multiple query files and databases.
    - Filters results based on the E-value threshold.
    - Outputs results in a structured CSV file with customizable column names.

Instructions:
    1. Input the path to the directory containing your BLAST databases.
    2. Choose the BLAST search method (e.g., blastn, tblastn).
    3. Input the path to the directory containing your query files.
    4. Specify where you want to store the results.
    5. Enter the maximum E-value for filtering results.

Example:
    Enter path to BLAST DBs: /path/to/blast_dbs
    What method would you like to use? (blastn, tblastn): blastn
    Enter path to query file(s): /path/to/query_files
    Where do you want to store the results?: /path/to/results
    Enter maximum E-value: 0.01

Author:
    Elijah R. Bring Horvath (https://github.com/ERBringHorvath)

License:
    This script is shared under MIT License, which allows for modification and redistribution with attribution.

Note:
    This script is intended for research and academic purposes. 
    Please ensure you have the necessary permissions to use the databases and query files with BLAST.

"""

