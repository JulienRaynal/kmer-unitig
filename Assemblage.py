import os
import sys
import argparse
import statistics
from Bio import SeqIO
from typing import Dict, List, Tuple

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


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
    #od = dict(sorted(kmer_dict.items()))
    print(f"Number of kmer made by jellyfish: {len(kmer_dict)}")
    return kmer_dict


def alignSequences(ordered_dict: Dict[str, int], kmer_size: int, unitig_size: int) -> List[Tuple[str,int]]:
    """
    Calls the search algorithm on the right and then the left side of the given kmer

    Cette fonction et celle en dessous peuvent être optimises mais par respect du temps et des regles je ne le ferais pas sur ce rendu
    Parameters
    ----------
    kmer: str
        A kmer
    ordered_dict: Dict[str, int]
        A dict of all the kmers as keys and their occurences as values
    """
    print("Building unitigs ... Please wait")
    unitigs: List[Tuple[str, int]] = []  # A list containing a dict with as key the full unitig and as value the average kmer count along the path
    for key, value in ordered_dict.copy().items():
        unitig: str = key
        support = value
        keys_to_delete: List[str] = []
        unitigR = buildUnitig(unitig, ordered_dict, kmer_size)
        unitigL = buildUnitig(unitig, ordered_dict, kmer_size, right_search=False)
        for i in range(0, len(unitigR), 1):
            unitig = unitig + unitigR[i][-1]
            if ordered_dict.get(unitigR[i]):
                existing_key = unitigR[i]
            else:
                existing_key = computeReverseComplement(unitigR[i])
            support = (support + ordered_dict.get(existing_key))/2
            keys_to_delete.append(existing_key)
        for i in range(len(unitigL) - 1, -1, -1):
            unitig = unitigL[i][0] + unitig
            if ordered_dict.get(unitigL[i]):
                existing_key = unitigL[i]
            else:
                existing_key = computeReverseComplement(unitigL[i])
            support = (support + ordered_dict.get(existing_key)) / 2
            keys_to_delete.append(existing_key)

        if len(unitig) >= unitig_size:
            unitigs.append((unitig, support))
            for key in keys_to_delete:
                del ordered_dict[key]
        else:
            continue
    return unitigs


def buildUnitig(unitig: str, ordered_dict: Dict[str, int], kmer_size: int, right_search=True) -> List[str]:
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
    unitig_list = [unitig]
    while not stop:
        complementary_chars: List[str] = []
        kmer_counts: List[int] = []

        if right_search is True:
            # get the prefix
            substring: str = unitig_list[-1][1:]
        else:  # get the suffix
            substring: str = unitig_list[0][:-1]

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
            if right_search:
                unitig_list.append(substring + complementary_chars[0])
            else:
                unitig_list.insert(0, complementary_chars[0] + substring)
            avg_support = (avg_support + kmer_counts[0] / 2)

    if right_search:
        unitig_list = unitig_list[1:]
    else:
        unitig_list = unitig_list[:-1]
    return unitig_list

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

def unitigs_to_fasta(unitigs: List[Tuple[str, int]]) -> None:
    """
    Write the unitigs as fasta
    Parameters
    ----------
    unitigs: List[Tuple[str, int]]
        A list containing tuples of sequences and average support

    Returns
    -------
    None
    """
    records: List[SeqRecord] = []
    for idx, unitig in enumerate(unitigs):
        record = SeqRecord(
            Seq(unitig[0]),
            id=f"unitig_{idx}_{len(unitig[0])}_{unitig[1]}",
            description=""
        )
        records.append(record)

    SeqIO.write(records, "unitigs.fasta", "fasta")

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
    if (args.inputfilename.split(".")[-1] == "fq" or args.inputfilename.split(".")[-1] == "fastq"):
        reads: list = cleanReads(reads, args.outputfilename, args.minreadquality, args.kmersize)

    #############################
    #       Creating kmer       #
    #############################
    kmer_dict: Dict[str, int] = jellyfish_kmer(args.outputfilename, args.kmersize)

    #############################
    #       Building unitig     #
    #############################
    unitigs = alignSequences(kmer_dict, args.kmersize, args.unitigsize)
    unitigs_to_fasta(unitigs)
