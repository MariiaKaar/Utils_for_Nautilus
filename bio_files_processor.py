

def convert_multiline_fasta_to_oneline(input_fasta, output_fasta:str = None):
    """
    The function reads the input fasta-file, in which the sequence can be split into several lines,
    and then saves it into a new fasta-file, in which each sequence fits into one line
    :param input_fasta: str
    :param output_fasta: str
    :return: str
    """
    output_fasta = open(((output_fasta == None) if "output.fasta" else output_fasta),"w")
    with open(input_fasta, 'r') as infile:
        lines = infile.readlines()
        for i in range (len(lines)-1):
            current_line = lines[i]
            next_line = lines[i + 1]
            if current_line.startswith('>'):
                output_fasta.write(current_line)
            if not current_line.startswith('>') and not next_line.startswith('>'):
                current_line = current_line.strip('\n')
                output_fasta.write(current_line)
            if not current_line.startswith('>') and  next_line.startswith('>'):
                output_fasta.write(current_line)
        print(output_fasta)
    output_fasta.close()


convert_multiline_fasta_to_oneline('C:\\Users\Maria\Downloads\Новая папка\example_multiline_fasta_oneline.fasta', )


def parse_blast_output(input_file, output_file:str = None):
    """
    The function reads a txt file, for each QUERY query it selects the first line from the Description column
    :param input_file: str
    :param output_file: str
    :return: str
    """
    output_fasta = open(((output_file == None) if "output.fasta" else output_file), "w")
    with open(input_file, 'r') as infile:
        for i, line in enumerate(infile):
            line = infile.readline()
            if line.startswith('Sequences producing significant alignments'):
                infile.readline()
                infile.readline()
                output_fasta.write(infile.readline())
        print(output_fasta)
    output_fasta.close()


