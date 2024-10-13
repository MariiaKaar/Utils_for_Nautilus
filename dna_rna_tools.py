#from modules import dna_rna_tools_modules
from modules import filter_fastq_modules


def run_dna_rna_tools(*args):
   pass

def filter_fastq(seq :dict[str, tuple[str,str]], gc_bounds:tuple=(0, 100), length_bounds:int=(0, 2 ** 32), quality_threshold: int=0 ):
    """
    Функция проверяет, что подаваемый словарь, состоящий из fastq-сиквенсов  соответствует условиям:
     1)находится в указанном диапазоне GC состава( если не указано, то все риды сохраняем)
     2)длина рида попадает в указанный интервал длины
     3)качество рида выше порогового

    :param fastq_file:dict
    :param gc_bounds:tuple
    :param length_bounds:
    :param quality_threshold:int
    :return:dict
    """
    filtered_fastq = {}
    for seq_name, seq_data in seq.items():
        seq, quality = seq_data

        gc_content = filter_fastq_modules.calculate_gc_bounds(seq) * 100
        seq_length = filter_fastq_modules.length_bounds(seq)
        q_score = filter_fastq_modules.quality_threshold(quality)

        if isinstance(gc_bounds, int):
            lower_gc_threshold, upper_gc_threshold = 0, gc_bounds
        else:
            lower_gc_threshold, upper_gc_threshold = gc_bounds

        if isinstance(length_bounds, int):
            lower_length_threshold, upper_length_threshold = 0, length_bounds
        else:
            lower_length_threshold, upper_length_threshold = length_bounds

        if (lower_gc_threshold <= gc_content <= upper_gc_threshold and
                lower_length_threshold <= seq_length <= upper_length_threshold and
                q_score >= quality_threshold):
            filtered_fastq[seq_name] = (seq, quality)

    return filtered_fastq

