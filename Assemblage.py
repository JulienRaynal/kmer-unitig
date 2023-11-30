import os
import sys
import argparse
import statistics
from Bio import SeqIO
from typing import Dict, List


def getFileExtension(path: str) -> str:
    """
    Get the extension of a file

    Parameters
    ----------
    path: str
        A path to a file
    
    Returns
    -------
    str
    """
    formats = {
        "fasta": "fasta",
        "fa": "fasta",
        "fastq": "fastq",
        "fq": "fastq"
    }

    file_extension = path.rsplit('.', 1)[-1]
    return formats[file_extension]


def getReads(path: str) -> iter:
    """
    Get all the reads from a file
    Manages only fasta or fastq

    Parameters
    ----------
    path: str
        A string to the path of the file to be read
    """

    # Get the extension of the file
    file_format = getFileExtension(path)
    return SeqIO.parse(path, file_format)


def cleanReads(records: iter, outputfile_name: str, min_quality: int = 30, kmer_size: int = 31) -> list:
    """
    Remove the reads by a minimal quality score

    Parameters
    ----------
    records: iter
        A iterable set of reads
    outputfile_name: str
        The name of the output file to save the reads
    min_quality: int (default = 0)
        The minimum average quality of a read to be kept
    kmer_size: The size of the wanted kmers

    Returns
    -----------
    List of kept reads
    """
    good_reads: list = []
    read_count: int = 0
    for record in records:
        if statistics.mean(record.letter_annotations["phred_quality"]) >= 30:
            good_reads.append(record)
            read_count += 1  # Doing this as we can't get the length of an iterable

    print(f"Found {len(good_reads)} reads with a quality above {min_quality}")

    getStats(good_reads, kmer_size)

    # Saving reads
    print("Saving", read_count, "reads")
    SeqIO.write(good_reads, outputfile_name, getFileExtension(outputfile_name))

    return good_reads


def getStats(records: list, kmer_size: int) -> None:
    """
    Prints the stats of the current list of reads

    Parameters
    ----------
    records: list
        A list of reads
    kmer_size: int
        The size of a kmer
    """
    total_length: int = 0

    for record in records:
        total_length += len(record.seq)

    print(f"The total length of all the reads is: {total_length}\n\
      The total number of reads is: {len(records)}\n\
      The total number of kmer is: {total_length - kmer_size + 1}\n\
      Theoritical upper bounds of kmers in the reads is: {total_length - len(records) * (kmer_size - 1)}")

    print(f"\nThe coverage is : {total_length / 10777}\n\
          The real number of kmer is : {10777 - kmer_size + 1}")


def jellyfish_kmer(output_file: str, kmer_size: int) -> Dict[str, int]:
    """
    Uses jellyfish to get the kmer of the reads

    Parameters
    ----------
    output_file: str
        The file containing the cleaned reads and also the file to save the jellyfish output
    kmer_size: int
        The wanted size for the kmers

    returns
    -------
    Ordered dict
    """
    print("Using jellyfish to create kmers")
    os.system(f"jellyfish count -m {kmer_size} -s 100M -t 4 -C {output_file} -o counted_kmers.jf")
    os.system(f"jellyfish dump counted_kmers.jf > {output_file} && rm counted_kmers.jf")

    kmer_dict: Dict[str, int] = {}
    print("Creating the kmer dict")
    for record in SeqIO.parse(output_file, getFileExtension(output_file)):
        kmer_dict[str(record.seq)] = int(record.id)
    od = dict(sorted(kmer_dict.items()))
    print(f"Number of kmer made by jellyfish: {len(od)}")
    return od


def alignSequence(kmer: str, ordered_dict: Dict[str, int], kmer_size):
    """
    Calls the search algorithm on the right and then the left side of the given kmer
    
    Parameters
    ----------
    kmer: str
        A kmer
    ordered_dict: Dict[str, int]
        A dict of all the kmers as keys and their occurences as values
    """
    unitig: str = kmer
    possibilities: List[str] = 'A', 'T', 'C', 'G'
    avg_support: int = 0
    stop: bool = False

    print(f"Searching the kmer aligning with {kmer} on the right")
    unitig = buildUnitig(unitig, ordered_dict, kmer_size)
    buildUnitig(unitig, ordered_dict, kmer_size, right_search=False)


def buildUnitig(unitig: str, ordered_dict: Dict[str, int], kmer_size: int, right_search=True) -> str:
    """
    Builds a unitig based on the kmer found in the dict built by jellifish
    Parameters
    ----------
    unitig
        A kmer to begin the unitig
    ordered_dict
        A dict holding all the kmer and the number of times it has appeared in jellyfish
    kmer_size
        The size of the kmer
    right_search
        If it tries to build the right side of the unitig

    Returns
    -------
    str
    """
    stop: bool = False
    avg_support: float = 0

    while not stop:
        complementary_chars: List[str] = []
        kmer_counts: List[int] = []

        if right_search is True:  # get the prefix
            substring: str = unitig[len(unitig) - kmer_size + 1:]
        else:  # get the suffix
            substring: str = unitig[:kmer_size - 1]

        reverse_complement_substring: str = computeReverseComplement(substring)
        # Check if another kmer can be added depending on the 4 possibilities
        for char in ['A', 'T', 'C', 'G']:
            kmer_to_search: str = substring + char if right_search else char + substring
            reverse_to_search: str = char + reverse_complement_substring if right_search else reverse_complement_substring + char
            findComplementaryKmer(kmer_to_search, ordered_dict, char, complementary_chars, kmer_counts)
            findComplementaryKmer(reverse_to_search, ordered_dict, char, complementary_chars, kmer_counts, reverse=True)
        # If there's more than one option we've reached a crossroad and we stop
        if len(complementary_chars) != 1:
            stop = True
        else:
            unitig = unitig + complementary_chars[0] if right_search else complementary_chars[0] + unitig
            print(f"Added a char to unitig: {unitig}")
            avg_support = (avg_support + kmer_counts[0] / 2)

    print(unitig)
    return unitig

def findComplementaryKmer(kmer_to_search: str, ordered_dict: Dict[str, int], char: str,
                          complementary_chars: List[str], kmer_counts: List[int], reverse=False):
    """
    Manages finding in a dict the occurences of a given kmer and add the data to the right list
    """
    kmer_count: int = ordered_dict.get(kmer_to_search)
    if kmer_count is not None:
        kmer_counts.append(kmer_count)
        if reverse:
            complementary_chars.append(computeReverseComplement(char))
        else:
            complementary_chars.append(char)


def computeReverseComplement(kmer: str) -> str:
    """
    Computes the reverse complement of the given sequence

    Parameters
    ----------
    kmer: str
        A kmer as string for which the complement needs to be computed

    Returns
    -------
    str
    """
    base_complement = {"A": "T",
                       "T": "A",
                       "C": "G",
                       "G": "C"}
    reverse: str = kmer[::-1]
    complement: str = ""
    for char in reverse:
        complement += base_complement[char]

    return complement


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Kmer Grapher")
    parser.add_argument("--inputfilename", help="The path of the input file", required=True, type=str)
    parser.add_argument("--outputfilename", help="The path for the output file", required=True, type=str)
    parser.add_argument("--kmersize", help="The size of kmers", choices=range(2, 32), required=True, type=int)
    parser.add_argument("--minkmercount", help="Minimum count of a kmer to be kept", required=True, type=int)
    parser.add_argument("--minunitiglength", help="Minimum length of the unitig", required=True, type=int)
    parser.add_argument("--minreadquality", help="The minimum quality a read should have on average", type=int,
                        default=30)
    parser.add_argument("--unitigsize", help="The size of the unitig", type=int, required=True)

    args = parser.parse_args()

    #############################
    #       Sorting reads       #
    #############################
    reads: iter = getReads(args.inputfilename)
    reads: list = cleanReads(reads, args.outputfilename, args.minreadquality, args.kmersize)

    #############################
    #       Creating kmer       #
    #############################
    ordered_dict: Dict[str, int] = jellyfish_kmer(args.outputfilename, args.kmersize)

    #############################
    #       Building unitig     #
    #############################
    unitigs: List[Dict[
        str, int]] = []  # A list containing a dict with as key the full unitig and as value the average kmer count along the path
    for i in range(len(list(ordered_dict.keys()))):
        alignSequence(list(ordered_dict.keys())[i], ordered_dict, args.kmersize)

    # d = dict = {
    # "AAAAAAATCGTTTCGGGATGATGCATAGCAT": 1,
    # "AAAAAATCGTTTCGGGATGATGCATAGCATA": 1,
    # "AAAAATCGTTTCGGGATGATGCATAGCATAT": 1,
    #        }
    # alignSequence(list(d.keys())[0], d, args.kmersize)
