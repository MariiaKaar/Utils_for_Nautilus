# Utils_for_Nautilus

**Utils_for_Nautilus** is a package for processing and filtering FASTQ files.

**Authors:**

* **Development**: Maria Kravchenko (https://github.com/MariiaKaar)
* **Supervisor, Concept**: Institute of Bioinformatics (https://bioinf.me/en),

  Saint Petersburg, Russia

## Content
* [Installation](#installation)
* [Usage Instructions](#usage-instructions)
* [Contacts](#contacts)
* [Functions](#functions)

## Usage Instructions

To work with the package, you need to import the main script and call any of the functions.

## Functions

The program contains two main functions and additional modules for them.

1. **run_dna_rna_tools** (this function will be implemented in the future)

2. **filter_fastq** accepts a dictionary with reads and quality metrics, and calculates:
   - the percentage of GC content in the read.
   - the threshold quality value of the read.
   - the length of the read.
   -  reading a file at a given path and translating it into a dictionary
   -  writes the sequences dictionary to the output_fastq file, and saves the data to the filtered folder
3. **bio_file_processor**
  - reads the input fasta-file, in which the sequence can be split into several lines,
    and then saves it into a new fasta-file, in which each sequence fits into one line
  - reads a txt file, for each QUERY query it selects the first line from the Description column


## Installation

To run Utils_for_Nautilus, download the folder to your computer or clone the repository.

```bash
git clone git@github.com:MariiaKaar/Utils_for_Nautilus.git

Contacts
My email: mariakr55@mail.ru
