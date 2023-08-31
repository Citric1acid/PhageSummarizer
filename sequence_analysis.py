# sequence analysis
# define some helping functions

import logging
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO


class PhageSequence:
    contig: str
    start: int
    end: int
    length: int
    info: str

    def __init__(self, contig, start, end, info=''):
        self.contig = contig
        self.start = start
        self.end = end
        self.length = end - start + 1
        self.info = info

    def __str__(self):
        return f'contig: {self.contig}, start: {self.start}, end: {self.end}, length: {self.length}, info: {self.info}'

    def __repr__(self):
        return f'{self.__class__.__name__}(contig={self.contig}, start={self.start}, end={self.end}, length={self.length}, info={self.info})'


def readfile(file, mode='r') -> str:
    with open(file, mode) as f:
        return f.read()


def read_sequences(file, format='fasta') -> list[SeqRecord]:
    seq = list(SeqIO.parse(file, format))
    logging.info(f'Reading {file} with {len(seq)} records')
    return seq


def contig_length(seqs: list[SeqRecord]) -> dict:
    return {seq.id: len(seq) for seq in seqs}


def clear_name(s: str) -> str:
    # remove path and suffix
    return s.split('/')[-1].removesuffix(".fasta")


def change_suffix(s: str) -> str:
    # change .fa or .fna into .fasta
    if s.endswith('.fa') or s.endswith('.fna'):
        return s.rsplit('.', 1)[0] + '.fasta'
    else:
        return s


def seqs_to_str(seqs: list[SeqRecord], format='fasta') -> str:
    handle = StringIO()
    SeqIO.write(seqs, handle, format)
    return handle.getvalue()


def count(lst, func):
    return len([x for x in lst if func(x)])


def find_item(lst, func):
    for x in lst:
        if func(x):
            return x


def to_table(lst) -> pd.DataFrame:
    return pd.DataFrame([obj.__dict__ for obj in lst])


def union_of_intervals(intervals) -> list:
    result = []
    for x in sorted(intervals):
        # put the smallest interval
        if not result:
            result.append(x)
            continue
        # add other intervals
        current = result[-1][1]
        if x[0] > current:
            # no intersection
            result.append(x)
        else:
            # x[0] is in the previous interval
            result[-1][1] = max(current, x[1])
    return result
