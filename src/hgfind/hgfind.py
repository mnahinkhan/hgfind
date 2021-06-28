#!/usr/bin/env python3
"""
Module for converting a gene name to its coordinates on hg38

Author: Nahin Khan <mnahinkhan@gmail.com>
"""

import argparse
import os
import pickle
import sys
from pathlib import Path
from typing import Dict, Union

PATH_TO_BIO_MART = Path(__file__).parent / "biomart-gene-coordinates.txt"
PATH_TO_PICKLE = Path(__file__).parent / "biomart-gene-coordinates.pickle"


class WrongGeneName(Exception):
    """
    Exception raised when a gene name is not recognized

    """


class Chromosome:
    """
    This class allows one to specify chromsomes from the human genome.
    Benefits include comparison of chromosomes that are sex chromosomes,
    autosomes, or genes included in the mitochondrial DNA.

    """

    def __init__(self, n: Union[int, str]):

        if isinstance(n, int) and 1 <= n <= 22:
            self.chr_n = n
        elif isinstance(n, str) and n.upper() in ("X", "Y", "MT", "M"):
            if n.upper() == "X":
                self.chr_n = 23
            elif n.upper() == "Y":
                self.chr_n = 24
            elif n.upper() == "MT":
                self.chr_n = 25
            elif n.upper() == "M":
                self.chr_n = 25
        elif isinstance(n, str) and is_int(n) and 1 <= int(n) <= 22:
            self.chr_n = int(n)
        else:
            raise ValueError("expected X, Y, M(T), or [1-22]")

    def __str__(self):
        return (
            str(self.chr_n)
            if self.chr_n <= 22
            else "X"
            if self.chr_n == 23
            else "Y"
            if self.chr_n == 24
            else "M"
        )

    def __lt__(self, other):
        return self.chr_n < other.chr_n

    def __eq__(self, other):
        if isinstance(other, str):
            return str(self) == other
        if isinstance(other, int):
            return self.chr_n == other
        if isinstance(other, Chromosome):
            return self.chr_n == other.chr_n
        raise TypeError

    def __hash__(self):
        return hash(self.chr_n)

    def __gt__(self, other):
        return self.chr_n > other.chr_n

    def __repr__(self):
        return self.__str__()

    def __int__(self):
        return self.chr_n


def hgfind(gene: str) -> Dict:
    """
    Given a string containing a gene name from the human genome, returns its
    location on hg38.

    :param gene: a string representing the name of a human gene.
    :returns: On success, a dictionary containing the following keys and
      associated values:
        'chr_n': the Chromosome object on which the gene lies
        'start_coord': the start coordinate of the gene on the chromosome
        'end_coord': the end coordinate of the gene on the chromosome
        'official_name': the standardized name for the specified gene
      On failure, an exception is raised

    """

    is_not_pickled = True
    if os.path.isfile(PATH_TO_PICKLE):
        is_not_pickled = False
        with open(PATH_TO_PICKLE, "rb") as handle:
            try:
                name_to_official, official_to_coord = pickle.load(handle)
            except ModuleNotFoundError:
                is_not_pickled = True
            except AttributeError:
                is_not_pickled = True
            except ValueError:
                is_not_pickled = True

    if is_not_pickled:
        name_to_official, official_to_coord = file_to_dicts(PATH_TO_BIO_MART)
        # Create parent dir if needed
        Path(PATH_TO_PICKLE).parent.mkdir(parents=True, exist_ok=True)
        with open(PATH_TO_PICKLE, "wb") as handle:
            pickle.dump(
                [name_to_official, official_to_coord],
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )

    gene = gene.upper()

    if (
        gene in name_to_official
        and name_to_official[gene] in official_to_coord
    ):
        rna_info = {}
        (
            rna_info["chr_n"],
            rna_info["start_coord"],
            rna_info["end_coord"],
        ) = official_to_coord[name_to_official[gene]]
        rna_info["official_name"] = name_to_official[gene]

    else:
        raise WrongGeneName(
            {"message": "The input gene could not be recognized", "gene": gene}
        )

    return rna_info


def file_to_dicts(path):
    """
    This function converts the BioMart file given by the path into two
    dictionaries: the first one takes a gene name as a key and gives as value
    the official symbol for the gene. The second one then takes the official
    symbol as key and gives back the chromosome number, the start coordinate,
    and the end coordinate as its value.

    :param path: path to BioMart file specifying gene coordinates
    :returns: two dictionaries as tuple
        (gene name -> official name, official name -> genome location)

    """

    name_to_official: Dict[str, str] = {}
    official_to_coord = {}

    with open(path) as bio_mart_file:
        # This first line just says "Gene start (bp), Gene end (bp), Gene name,
        # Gene Synonym, WikiGene name, UniProtKB Gene Name symbol,
        # Chromosome/Scaffold Name"
        bio_mart_file.readline()

        # To get to an actual line:
        line = bio_mart_file.readline()
        while line:

            (
                start_coord,
                end_coord,
                gene_name,
                gene_syn_name,
                wiki_name,
                uniprot_name,
                str_chr_no,
            ) = [wss.strip() for wss in line.split("\t")]

            names = [
                name
                for name in [gene_name, gene_syn_name, wiki_name, uniprot_name]
                if name != ""
            ]

            if str_chr_no in ("X", "Y", "MT"):
                chr_no = Chromosome(str_chr_no)

            elif len(str_chr_no) <= 2:
                chr_no = Chromosome(int(str_chr_no))

            elif "HSCHR" in str_chr_no.upper():
                str_chr_no = str_chr_no[
                    str_chr_no.upper()
                    .find("HSCHR") : str_chr_no.upper()
                    .find("HSCHR")
                    + 7
                ]

                str_chr_no = str_chr_no[5:]

                if str_chr_no[1].isdigit():
                    chr_no = Chromosome(int(str_chr_no))

                elif str_chr_no[0] == "X" or str_chr_no[0] == "Y":
                    chr_no = Chromosome(str_chr_no[0])

                elif str_chr_no == "MT":
                    chr_no = Chromosome(str_chr_no)

                else:
                    chr_no = Chromosome(int(str_chr_no[0]))

            else:
                line = bio_mart_file.readline()
                continue

            start_coord = int(start_coord)
            end_coord = int(end_coord)

            # Maybe one of the names is already in the dictionary.Then, link all of
            # the names to the official title for that name. Otherwise, let
            # "gene name" be official.
            official_name = name_to_official.get(gene_name, gene_name)

            for name in names:
                if name in name_to_official:
                    official_name = name_to_official[name]
                    if official_name not in official_to_coord:
                        continue

                    (
                        prev_chr_no,
                        prev_start_coord,
                        prev_end_coord,
                    ) = official_to_coord[official_name]

                    if prev_chr_no != chr_no:
                        name_to_official[name] = name
                        official_name = name
                    continue

                name_to_official[name] = official_name

            # Now get the start and end coord of this row:
            (
                prev_chr_no,
                prev_start_coord,
                prev_end_coord,
            ) = official_to_coord.get(official_name, (-1, [], []))

            if int(prev_chr_no) != -1 and prev_chr_no != chr_no:
                line = bio_mart_file.readline()
                continue

            new_chr_no = chr_no
            new_start_coord = prev_start_coord + [start_coord]
            new_end_coord = prev_end_coord + [end_coord]

            official_to_coord[official_name] = (
                new_chr_no,
                new_start_coord,
                new_end_coord,
            )

            line = bio_mart_file.readline()

    print(
        "Generating pickle (cache) file for faster lookup next time...",
        file=sys.stderr,
    )
    for official_name, coords in official_to_coord.items():
        chr_no, start_coord_array, end_coord_array = coords

        start_coord = min(start_coord_array)
        end_coord = max(end_coord_array)

        official_to_coord[official_name] = (chr_no, start_coord, end_coord)

    return name_to_official, official_to_coord


def is_int(in_obj):
    """
    Checks if the input represents an integer and returns true iff so

    """
    try:
        int(in_obj)
        return True
    except ValueError:
        return False


def main():
    """
    Parses command line gene argument and prints result on stdout

    """
    parser = argparse.ArgumentParser(
        description="Get the human genome (hg38) coordinates of a gene",
    )
    parser.add_argument("gene", help="The name of a gene")
    args = parser.parse_args()

    test_gene = args.gene
    try:
        result = hgfind(test_gene)
    except WrongGeneName:
        print(f"{test_gene} not recognized as a gene")
        sys.exit(1)

    test_gene_name = result["official_name"]
    chr_n = result["chr_n"]
    start = result["start_coord"]
    end = result["end_coord"]
    print(f"{test_gene_name} => {chr_n}:{start}-{end}")


if __name__ == "__main__":
    main()
