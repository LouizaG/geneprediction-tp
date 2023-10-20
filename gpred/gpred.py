import argparse
import sys
import os
import csv
import re
import textwrap
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage="{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument(
        '-g',
        dest='min_gene_len',
        type=int,
        default=50,
        help="Minimum gene length to consider (default 50).")
    parser.add_argument(
        '-s',
        dest='max_shine_dalgarno_distance',
        type=int,
        default=16,
        help="Maximum distance from start codon "
        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument(
        '-d',
        dest='min_gap',
        type=int,
        default=40,
        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path,
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome.
    """
    with open(fasta_file, "r") as file:
        return ''.join(line.strip().upper()
                       for line in file if not line.startswith(">"))


def find_start(start_regex: Pattern, sequence: str,
               start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None.
    """
    match = start_regex.search(sequence, start, stop)
    return match.start(0) if match else None


def find_stop(stop_regex: Pattern, sequence: str,
              start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None.
    """
    match_object = stop_regex.finditer(sequence, start)
    for match in match_object:
        if match.start(0) % 3 == start % 3:
            position = match.start(0)
            return position
    return None


def has_shine_dalgarno(
        shine_regex: Pattern,
        sequence: str,
        start: int,
        max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    window_start = start - max_shine_dalgarno_distance
    window_end = start - 6

    # If start of the window is negative, return False
    if window_start < 0:
        return False
    # Search for the Shine-Dalgarno sequence within the window
    match = shine_regex.search(sequence, window_start, window_end)

    return bool(match)


def predict_genes(sequence: str,
                  start_regex: Pattern,
                  stop_regex: Pattern,
                  shine_regex: Pattern,
                  min_gene_len: int,
                  max_shine_dalgarno_distance: int,
                  min_gap: int) -> List[List[int]]:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Maximum distance of the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    position_courante = 0
    predicted_genes = []

    while len(sequence) - position_courante >= min_gap:
        start_pos = find_start(
            start_regex,
            sequence,
            position_courante,
            len(sequence))
        if start_pos is not None:
            stop_pos = find_stop(stop_regex, sequence, start_pos)
            if stop_pos is not None and stop_pos - start_pos >= min_gene_len:
                if has_shine_dalgarno(
                        shine_regex,
                        sequence,
                        start_pos,
                        max_shine_dalgarno_distance):
                    predicted_genes.append([start_pos + 1, stop_pos + 3])
                    position_courante = stop_pos + 3 + min_gap
                    continue
        position_courante += 1

    return predicted_genes


def write_genes_pos(predicted_genes_file: Path,
                    probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def write_genes(fasta_file: Path,
                sequence: str,
                probable_genes: List[List[int]],
                sequence_rc: str,
                probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i, gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i + 1, os.linesep,
                    textwrap.fill(sequence[gene_pos[0] - 1:gene_pos[1]])))
            i = i + 1
            for j, gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i + 1 + j, os.linesep,
                            textwrap.fill(sequence_rc[gene_pos[0] - 1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    # AGGA ou GGAGG
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')
    # Arguments
    args = get_arguments()
    sequence = read_fasta(args.genome_file)
    find_start_result_updated = find_start(
        start_regex, sequence, 0, len(sequence))
    print(f"codon start: ", find_start_result_updated)
    find_stop_result = find_stop(
        stop_regex, sequence, find_start_result_updated)
    print(f"codon stop: ", find_stop_result)
    has_shine_dalgarno_test = has_shine_dalgarno(
        shine_regex,
        sequence,
        find_start_result_updated,
        args.max_shine_dalgarno_distance)
    print("has shine dalgarno ?:", has_shine_dalgarno_test)

    predicted_genes_5_to_3 = predict_genes(
        sequence,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap)

    # We reverse and complement
    reversed_sequence = reverse_complement(sequence)
    # Predict genes in reversed sequence (3' to 5' direction)

    predicted_genes_3_to_5_reversed = predict_genes(
        reversed_sequence,
        start_regex,
        stop_regex,
        shine_regex,
        args.min_gene_len,
        args.max_shine_dalgarno_distance,
        args.min_gap)

    # Correct the positions for 3' to 5' predictions to make them appear in 5'
    # to 3' direction

    predicted_genes_3_to_5 = []
    for gene in predicted_genes_3_to_5_reversed:
        start, end = gene
        corrected_start = len(sequence) - end + 1
        corrected_end = len(sequence) - start + 1
        predicted_genes_3_to_5.append([corrected_start, corrected_end])

    # Combine results
    all_predicted_genes = predicted_genes_5_to_3 + predicted_genes_3_to_5

    # return sorted(all_predicted_genes, key=lambda x: x[0])
    write_genes_pos(args.predicted_genes_file, all_predicted_genes)
    write_genes(
        args.fasta_file,
        sequence,
        predicted_genes_5_to_3,
        reversed_sequence,
        predicted_genes_3_to_5_reversed)


if __name__ == '__main__':
    main()
