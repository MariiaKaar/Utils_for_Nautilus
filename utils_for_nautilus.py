from abc import ABC, abstractmethod
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio import SeqIO
import logging
import argparse
from typing import Union, Tuple


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

    def __init__(self, seq, dict):
        self.seq = seq
        self.dict = dict

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        return f"your sequence: {self}"

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
        super().__init__(seq, dna_complement_dict)

    def transcribe(self):
        return self.seq.replace("T", "U").replace("t", "u")


class RNASequence(NucleicAcidSequence):
    """
    Represents an RNA sequence.

    Provides methods specific to RNA sequences.
    """

    def __init__(self, seq):
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
        super().__init__(seq, rna_complement_dict)


class AminoAcidSequence(BiologicalSequence):
    """
    Represents an amino acid sequence.

    Provides methods for working with amino acid sequences, including calculating the molecular weight.
    """

    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)

    def __repr__(self):
        print(f"your sequence: {self}")

    def check_alphabet(self):
        valid_bases = set("ACDEFGHIKLMNPQRSTVWYBZX*")
        return all(base.upper() in valid_bases for base in self.seq)

    def __getitem__(self, index):
        return self.seq[index]

    def __str__(self):
        return self.seq

    def calculate_molecular_weight(self):
        return len(self.seq) * 110


logging.basicConfig(
    filename="filter_fastq.log",
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def filter_fastq(input_file, output_file, gc_bounds, length_bounds, quality_threshold):
    """ "Filters reads from a FASTQ file based on GC content, read length, and average quality score.

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
    try:
        if isinstance(gc_bounds, int):
            lower_gc_threshold, upper_gc_threshold = 0, gc_bounds
        else:
            lower_gc_threshold, upper_gc_threshold = gc_bounds

        if isinstance(length_bounds, int):
            lower_length_threshold, upper_length_threshold = 0, length_bounds
        else:
            lower_length_threshold, upper_length_threshold = length_bounds

        logging.info(f"Start filtering: {input_file} -> {output_file}")
        logging.info(
            f"GC bounds: {gc_bounds}, Length bounds: {length_bounds}, Quality threshold: {quality_threshold}"
        )

        good_reads = (
            rec
            for rec in SeqIO.parse(input_file, "fastq")
            if (
                lower_gc_threshold <= gc_fraction(rec.seq) * 100 <= upper_gc_threshold
                and lower_length_threshold <= len(rec.seq) <= upper_length_threshold
                and sum(rec.letter_annotations["phred_quality"]) / len(rec.seq)
                >= quality_threshold
            )
        )

        count = SeqIO.write(good_reads, output_file, "fastq")
        logging.info(f"Filtering end. Write {count} reads in {output_file}")
        return count
    except Exception as e:
        logging.error(f"Error in filtering: {str(e)}")
        raise




def parse_gc_bounds(value: str) -> Union[int, Tuple[int, int]]:
    """Parse argument - gc_bounds"""
    if "," in value:
        lower, upper = map(int, value.split(","))
        return (lower, upper)
    return int(value)


def parse_length_bounds(value: str) -> Union[int, Tuple[int, int]]:
    """Parse argument - length_bounds"""
    if "," in value:
        lower, upper = map(int, value.split(","))
        return (lower, upper)
    return int(value)


def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTQ files by GC content, length and quality",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("input_file", help="Input FASTQ file path")
    parser.add_argument("output_file", help="Output FASTQ file path")
    parser.add_argument(
        "--gc",
        dest="gc_bounds",
        type=parse_gc_bounds,
        required=True,
        help="GC bounds as integer  or range (",
    )
    parser.add_argument(
        "--length",
        dest="length_bounds",
        type=parse_length_bounds,
        required=True,
        help="Length bounds as integer  or range ",
    )
    parser.add_argument(
        "--quality", type=float, required=True, help="Minimum average quality threshold"
    )

    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()

    if args.verbose:
        print(f"Filtering {args.input_file} with parameters:")
        print(f"GC bounds: {args.gc_bounds}")
        print(f"Length bounds: {args.length_bounds}")
        print(f"Quality threshold: {args.quality}")

    count = filter_fastq(
        input_file=args.input_file,
        output_file=args.output_file,
        gc_bounds=args.gc_bounds,
        length_bounds=args.length_bounds,
        quality_threshold=args.quality,
    )

    if args.verbose:
        print(f"Saved {count} reads to {args.output_file}")


if __name__ == "__main__":
    try:
        filter_fastq(
            input_file="example_fastq.fastq",
            output_file="output.fastq",
            gc_bounds=(30, 70),
            length_bounds=(50, 150),
            quality_threshold=20.0,
        )
    except Exception as e:
        print(f"An error occurred: {e}")
