import pytest
from ..utils_for_nautilus import filter_fastq
import os
from Bio import SeqIO
from Bio.Seq import Seq


@pytest.fixture
def create_test_fastq(tmp_path):
    """Fixtura creates a test FASTQ file with different sequence characteristics"""
    test_fastq = tmp_path / "input.fastq"
    test_data = [
        (Seq("ACGT"), [40, 40, 40, 40], "read1"),  # 50% GC, high quality
        (Seq("AAAA"), [10, 10, 10, 10], "read2"),  # 0% GC, low quality
        (Seq("GGGG"), [30, 30, 30, 30], "read3"),  # 100% GC, average quality
        (Seq("ACGTACGT"), [25] * 8, "read4"),  # 50% GC,
    ]

    records = [
        SeqIO.SeqRecord(seq=seq, id=seq_id, letter_annotations={"phred_quality": qual})
        for seq, qual, seq_id in test_data
    ]
    SeqIO.write(records, test_fastq, "fastq")
    yield test_fastq


class TestBasicFunctionality:
    """Tests of basic functionality"""

    def test_file_creation(self, create_test_fastq, tmp_path):
        output_file = tmp_path / "output.fastq"
        count = filter_fastq(
            input_file=str(create_test_fastq),
            output_file=str(output_file),
            gc_bounds=100,
            length_bounds=10,
            quality_threshold=0,
        )
        assert os.path.exists(output_file)
        assert count == 4


class TestFiltering:
    """Filtering tests with parameterisation"""

    @pytest.mark.parametrize(
        "gc_bounds, expected_ids",
        [
            ((0, 100), ["read1", "read2", "read3", "read4"]),
            ((0, 0), ["read2"]),  # Only 0% GC quality
            ((100, 100), ["read3"]),  # GC quality 100%
            ((30, 70), ["read1", "read4"]),  # Medium GC-composition
            (50, ["read2", "read1", "read4"]),  # GC as int (0-50%)
        ],
    )
    def test_gc_filter(self, create_test_fastq, tmp_path, gc_bounds, expected_ids):
        output_file = tmp_path / "output.fastq"
        filter_fastq(
            input_file=str(create_test_fastq),
            output_file=str(output_file),
            gc_bounds=gc_bounds,
            length_bounds=10,
            quality_threshold=0,
        )

        records = list(SeqIO.parse(output_file, "fastq"))
        assert {rec.id for rec in records} == set(expected_ids)

    @pytest.mark.parametrize(
        "length_bounds, expected_ids",
        [
            ((1, 10), ["read1", "read2", "read3", "read4"]),
            ((4, 4), ["read1", "read2", "read3"]),
            ((5, 8), ["read4"]),
            ((0, 3), []),  # length_bounds 0-3
        ],
    )
    def test_length_filter(
        self, create_test_fastq, tmp_path, length_bounds, expected_ids
    ):
        output_file = tmp_path / "output.fastq"
        filter_fastq(
            input_file=str(create_test_fastq),
            output_file=str(output_file),
            gc_bounds=100,
            length_bounds=length_bounds,
            quality_threshold=0,
        )

        records = list(SeqIO.parse(output_file, "fastq"))
        assert {rec.id for rec in records} == set(expected_ids)


class TestErrorHandling:
    """Error handling tests"""

    @pytest.mark.parametrize(
        "invalid_params",
        [
            {"gc_bounds": (100, 0)},  # Reverse GC boundaries
            {"length_bounds": (10, 1)},  # Reverse length limits
            {"quality_threshold": -10},  # Negative quality
        ],
    )
    def test_invalid_parameters(self, create_test_fastq, tmp_path, invalid_params):
        output_file = tmp_path / "output.fastq"
        with pytest.raises(ValueError):
            filter_fastq(
                input_file=str(create_test_fastq),
                output_file=str(output_file),
                gc_bounds=invalid_params.get("gc_bounds", 100),
                length_bounds=invalid_params.get("length_bounds", 10),
                quality_threshold=invalid_params.get("quality_threshold", 0),
            )
