# GC counting.
def calculate_gc_bounds(seq):
    """
     Функция считает долю GC в риде
    :param seq: str
    :return: int,float
    """
    sum_g = seq.lower().count("g")
    sum_c = seq.lower().count("c")
    gc_content = (sum_g + sum_c) / length_bounds(seq)
    return gc_content


# length calculation.
def length_bounds(seq):
    """
    Функция считает длину рида
    :param seq: str
    :return: int
    """
    seq_length = len(seq)
    return seq_length


# Sequencing Quality Scores.
def quality_threshold(quality):
    """
    Функция считает пороговое значение среднего качества рида
    :param quality:int
    :return:int
    """
    seq_length = length_bounds(quality)
    summ_q_score = 0
    for letter in quality:
        letter_q_score = ord(letter) - 33
        summ_q_score += letter_q_score

    q_score = summ_q_score / seq_length
    return q_score

