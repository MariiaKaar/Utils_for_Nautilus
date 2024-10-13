import os
import shutil



def calculate_gc_bounds(seq):
    """
     The function counts the GC share in the sequencer
    :param seq: str
    :return: int,float
    """
    sum_g = seq.lower().count("g")
    sum_c = seq.lower().count("c")
    gc_content = (sum_g + sum_c) / length_bounds(seq)
    return gc_content



def length_bounds(seq):
    """
    The function counts the length of the sequences
    :param seq: str
    :return: int
    """
    seq_length = len(seq)
    return seq_length



def quality_threshold(quality):
    """
    The function counts the threshold value of the average quality of the sequences
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


def fastq_to_dict(input_fastq: str):
    """
    This function is for reading a file at a given path and translating it into a dictionary
    :param input_fastq: str
    :return:dict
    """
    seq = {}
    with open(input_fastq, 'r') as file:
        while True:
            # read 4 lines
            header = file.readline().strip()
            if not header:  # If the end of the file is reached
                break
            sequence = file.readline().strip()
            file.readline()  # Skip the comment line
            quality = file.readline()

            # save in dictionary
            seq[header] = (sequence, quality)
        print(seq)





def dict_to_fastq(seq: dict, input_fastq: str, output_fastq: str = None):
    """
    This function writes the sequences dictionary to the output_fastq file, and saves the data to the filtered folder
    :param seq: dict
    :param input_fastq: str
    :param output_fastq: str
    :return: str
    """
    output_fastq = open(((output_fastq == None) if "output.fastq" else output_fastq), "w")
    with open(input_fastq, 'r') as infile:  # Read the contents of a file without modifying the original file
        lines = infile.readline()
        for header, (sequence, quality) in seq.items():
            output_fastq.write(f"{header}\n")
            output_fastq.write(f"{sequence}\n")
            for i in range(2, len(lines), 3):
                output_fastq.write(lines[i])
            output_fastq.write(f"{quality}\n")
        print(output_fastq)
    output_fastq.close()
    path_to_cur_dir = os.getcwd()
    output_file_path = os.path.join(path_to_cur_dir,"filtered","output_fastq")
    shutil.move(path_to_cur_dir, output_file_path)





