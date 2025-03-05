# Utils_for_Nautilus
 **Utils_for_Nautilus** -are packages for processing and filtering FASTQ files
 
 **Authors:**
 
 * **Development** : Kravchenko Maria Antonovna( https://github.com/MariiaKaar)
 * **Leader, Idea**: Institute of Bioinformatics( https://bioinf.me/en),
    
   St. Petersburg, Russia
## Content
* [Installation](#Installation)
* [Instructions_for_Use](#Instructions_for_Use)
* [Classes](#Classes)
* [Functions](#Functions)
* [Contacts](#Contacts)

## 
 ## Installation

To run Utils_for_Nautilus, download the folder to your computer or clone the repository.

```bash
git clone git@github.com:MariiaKaar/Utils_for_Nautilus.git
```
## Instructions for use

To work with the package you must import the main script

## Classes

1. **BiologicalSequence** (Abstract Base Class)

   Provides a common interface for biological sequences  

   Key Methods:

  `__len__`: Returns the length of the sequence.

   `__repr__`: Returns a string representation of the sequence.

   check_alphabet: Validates if the sequence uses a valid alphabet.

   `__getitem__`: Allows indexing into the sequence.

   `__str__`: Returns the sequence as a string.

2. **NucleicAcidSequence**

   Key Methods:

   complement: Returns the complementary sequence.

   reverse: Returns the reversed sequence.

   reverse_complement: Returns the reverse complement of the sequence.

   check_alphabet: Validates if the sequence contains valid nucleic acid bases (A, T, G, C, U).

3. **DNASequence**

      Represents a DNA sequence.

      Key Methods:

      transcribe: Transcribes DNA to RNA
4. **RNASequence**

      Represents an RNA sequence

      Key Methods:

      Inherits all methods from NucleicAcidSequence

5. **AminoAcidSequence**

      Key Methods:

      check_alphabet: Validates if the sequence contains valid amino acid 
      
      symbols.

      calculate_molecular_weight: Estimates the molecular weight of the
      
      sequence.

## Functions


1. **filter_fastq** 

   Filters reads from a FASTQ file based on GC content, read length, and 

   average quality score.

   Inputs:

   input_file: Path to the input FASTQ file.

   output_file: Path to save the filtered FASTQ file.

   gc_bounds: Tuple or int for GC content bounds (in percentage).

   length_bounds: Tuple or int for read length bounds.

   quality_threshold: Minimum average quality score.


## Contacts

  My e-mail (mariakr55@mial.ru )
