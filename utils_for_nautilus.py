from abc import ABC, abstractmethod
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio import SeqIO

class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.

    Provides a common interface for different types of biological sequences.
    """
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def check_alphabet(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

class NucleicAcidSequence(BiologicalSequence):
    """
    Represents a nucleic acid sequence.

    Provides methods for working with nucleic acid sequences, including calculating the complement and reverse complement.
    """
    def __init__(self,seq,dict):
        self.seq = seq
        self.dict = dict
    
    def __len__(self):
        return len(self.seq)
    
    def __repr__(self):
        return f'your sequence: {self}' 

    def check_alphabet(self):
        valid_bases = set("ATGCUatgcu")
        return all(base in valid_bases for base in self.seq)

    def __getitem__(self, index):
        return self.seq[index]
    
    def __str__(self):
        return self.seq

    def complement(self):
        result = ""
        for letter in self.seq:
            result += self.dict[letter] 
        return result
    
    def reverse(self):
        return self.seq[::-1]
    
    def reverse_complement(self):
        complemented = self.complement()
        return complemented[::-1]
    
class DNASequence(NucleicAcidSequence):
    """
    Represents a DNA sequence.

    Provides methods specific to DNA sequences, including transcription.
    """
    def __init__(self, seq):
        dna_complement_dict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "a": "t",
        "t": "a",
        "g": "c",
        "c": "g",
    }
        super().__init__(seq,dna_complement_dict)
        
    def transcribe(self):
        return self.seq.replace('T', 'U').replace('t', 'u')
    
class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence.

    Provides methods specific to RNA sequences.
    """
    def __init__(self,seq):
        rna_complement_dict = {
        "A": "U",
        "U": "A",
        "G": "C",
        "C": "G",
        "a": "u",
        "u": "a",
        "g": "c",
        "c": "g",
    }
        super().__init__(seq,rna_complement_dict)


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence.

    Provides methods for working with amino acid sequences, including calculating the molecular weight.
    """
    def __init__(self,seq):
        self.seq = seq
    
    def __len__(self):
        return len(self.seq)
    
    def __repr__(self):
        print(f'your sequence: {self}')

    def check_alphabet(self):
        valid_bases = set("ACDEFGHIKLMNPQRSTVWYBZX*")
        return all(base.upper() in valid_bases for base in self.seq)

    def __getitem__(self, index):
        return self.seq[index]
    
    def __str__(self):
        return self.seq
    
    def calculate_molecular_weight(self):
        return len(self.seq) * 110 



def filter_fastq(input_file, output_file, gc_bounds, length_bounds , quality_threshold):
    """"Filters reads from a FASTQ file based on GC content, read length, and average quality score.

    Parameters:

    input_file : str
        Path to the input FASTQ file 

    output_file : str
        Path to the output FASTQ file where filtered reads will be saved.

    gc_bounds : tuple or int
        Bounds for filtering based on GC content (in percentage).
        - If an integer is provided, the bounds are treated as (0, gc_bounds).
        - If a tuple is provided, the bounds are treated as (lower_gc, upper_gc).

    length_bounds : tuple or int
        Bounds for filtering based on read length.
        - If an integer is provided, the bounds are treated as (0, length_bounds).
        - If a tuple is provided, the bounds are treated as (lower_length, upper_length).
       
    quality_threshold : float
        Threshold for the average Phred quality score.
        Reads with an average quality score below this threshold will be filtered out.
"""
    
    if isinstance(gc_bounds, int):
        lower_gc_threshold, upper_gc_threshold = 0, gc_bounds
    else:
        lower_gc_threshold, upper_gc_threshold = gc_bounds


    if isinstance(length_bounds, int):
        lower_length_threshold, upper_length_threshold = 0, length_bounds
    else:
        lower_length_threshold, upper_length_threshold = length_bounds


    good_reads = (
        rec
        for rec in SeqIO.parse(input_file, "fastq")
        if (lower_gc_threshold <= gc_fraction(rec.seq) * 100 <= upper_gc_threshold and
            lower_length_threshold <= len(rec.seq) <= upper_length_threshold and
            sum(rec.letter_annotations["phred_quality"]) / len(rec.seq) >= quality_threshold))


    count = SeqIO.write(good_reads, output_file, "fastq")
    return count