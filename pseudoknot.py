#! /usr/bin/python3
import glob
import itertools
import json
import logging
import os.path
import unittest
from dataclasses import dataclass
from typing import List


@dataclass
class Strand:
    '''
    A continuous fragment of RNA structure.

    When Strand object represents 5'-3' direction, then begin < end. Otherwise it is valid to have begin > end.

    Attributes:
        begin (int):        Index of the first nucleotide.
        end (int):          Index of the last nucleotide.
        sequence (str):     Strand sequence.
    '''
    begin: int
    end: int
    sequence: str

    def __eq__(self, other):
        return self.begin == other.begin and self.end == other.end and self.sequence == other.sequence

    def __str__(self):
        return '{:4d} {} {:4d}'.format(self.begin, self.sequence, self.end)


@dataclass
class Hairpin(Strand):
    '''
    A single strand of RNA structure with all nucleotides unpaired except the first and last which are paired together.

    In dot-bracket notation hairpins look like this (with variable number of dots):
    (...)
    '''


@dataclass
class Stem:
    '''
    Two strands of the same length, with all nucleotides forming pairs.

    In RNA structures, stems are anti-parallel. Therefore, the first strand is in 5'-3' direction, while the second
    strand in 3'-5'.

    In dot-bracket notation stems look like this (with variable number of brackets):
    ((((
    ))))

    Attributes:
        strand1 (Strand):   First instance of Strand.
        strand2 (Strand):   Second instance of Strand.
    '''
    strand1: Strand
    strand2: Strand

    def __eq__(self, other):
        return self.strand1 == other.strand1 and self.strand2 == other.strand2

    def __repr__(self):
        return '{} {} {}'.format(self.strand1.begin, self.strand2.begin,
                                 self.strand1.end - self.strand1.begin + 1)

    def __str__(self):
        return '{}\n     {}\n{}'.format(
            self.strand1, '|' * (self.strand1.end - self.strand1.begin + 1),
            self.strand2)

    def forms_pseudoknot_with(self, other):
        '''
        Check if this stem is in pseudoknot relation with another one.

        If stem1 is (i j len1) and stem2 is (k l len2), then these stems form a pseudoknot if i < k < j < l or
        k < i < l < j.

        Parameters:
            other (Stem):   Another instance of Stem to check against.

        Returns:
            bool:           Status of the pseudoknot relation.
        '''
        print("Dfsdfdsf")
        # TODO: implement this (change the False below to a meaningful expression)
        return False


@dataclass
class Pseudoknot:
    '''
    Two stems with pseudoknot relation between them.

    Attributes:
        stem1 (Stem):   First instance of Stem.
        stem2 (Stem):   Second instance of Stem.
    '''
    stem1: Stem
    stem2: Stem


class Encoder(json.JSONEncoder):
    '''
    Encoder for Strand, Hairpin, Stem and Pseudoknot objects to JSON.
    '''

    def default(self, o):
        if isinstance(o, Hairpin):
            return {
                'hairpin': {
                    'begin': o.begin,
                    'end': o.end,
                    'sequence': o.sequence
                }
            }
        if isinstance(o, Strand):
            return {
                'strand': {
                    'begin': o.begin,
                    'end': o.end,
                    'sequence': o.sequence
                }
            }
        if isinstance(o, Stem):
            return {
                'stem': {
                    'strand1': o.strand1,
                    'strand2': o.strand2
                }
            }
        if isinstance(o, Pseudoknot):
            return {
                'pseudoknot': {
                    'stem1': o.stem1,
                    'stem2': o.stem2
                }
            }


class Decoder(json.JSONDecoder):
    '''
    Decoder from JSON into Strand, Hairpin, Stem or Pseudoknot
    '''

    def __init__(self, *args, **kwargs):
        json.JSONDecoder.__init__(self, object_hook=self.object_hook, *args, **kwargs)

    def object_hook(self, dct):
        if 'strand' in dct:
            return dct['strand']
        if 'hairpin' in dct:
            strand = dct['hairpin']
            return Hairpin(strand.begin, strand.end, strand.sequence)
        if 'stem' in dct:
            return dct['stem']
        if 'pseudoknot' in dct:
            return dct['pseudoknot']
        if all(key in dct for key in ('begin', 'end', 'sequence')):
            return Strand(dct['begin'], dct['end'], dct['sequence'])
        if all(key in dct for key in ('strand1', 'strand2')):
            return Stem(dct['strand1'], dct['strand2'])
        if all(key in dct for key in ('stem1', 'stem2')):
            return Pseudoknot(dct['stem1'], dct['stem2'])


class BPSEQ:
    '''
    RNA structure.

    Attributes:
        entries (tuple):    A list of triples (i, s, j), where 'i' is the nucleotide index (starting from 1),
                            's' is a character in sequence and `j` is the index of paired nucleotide (or zero if
                            unpaired)
    '''

    @staticmethod
    def from_file(bpseq_path):
        '''
        Read a BPSEQ file and parse its content.

        Parameters:
            bpseq_path (str):   Path to a BPSEQ file.

        Returns:
            BPSEQ:              An instance of this class.
        '''
        entries = []
        with open(bpseq_path) as fd:
            for line in fd:
                fields = line.strip().split()
                if len(fields) != 3:
                    logging.warning(
                        'Failed to find 3 columns in BPSEQ line: {}', line)
                    continue
                entry = (int(fields[0]), fields[1], int(fields[2]))
                entries.append(entry)
        return BPSEQ(entries)

    def __init__(self, entries):
        self.entries = tuple((i, c, j) for i, c, j in entries)

    def __str__(self):
        return '\n'.join(
            ('{} {} {}'.format(i, c, j) for i, c, j in self.entries))

    def __eq__(self, other):
        return len(self.entries) == len(other.entries) and all(
            ei == ej for ei, ej in zip(self.entries, other.entries))

    def stems(self) -> List[Stem]:
        '''
        Find stem motifs.

        Returns:
            list:       A list of Stem objects.
        '''
        # TODO: implement this
        stems = []
        return stems

    def hairpins(self) -> List[Hairpin]:
        '''
        Find hairpin motifs.

        Returns:
            list:       A list of Hairpin objects.
        '''
        hairpins = []
        index = 0
        while index <= len(self.entries):
            if self.entries[index][2] == 0:
                index = index + 1
            else:
                sequence = ""
                begin = index + 1
                for i in range(index, len(self.entries)):
                    index = self.entries[index][0]
                    if(self.entries[index][2] == 0):
                        end = self.entries[index][0] - 1
                        break
                    else:
                        sequence = sequence + self.entries[index][1]
                hairpins.append(Hairpin(begin, end, sequence))
        return hairpins

    def pseudoknots(self) -> List[Pseudoknot]:
        '''
        Find pseudoknots.

        Returns:
            list:       A list of Pseudoknot objects.
        '''
        pseudoknots = []
        for stem1, stem2 in itertools.combinations(self.stems(), 2):
            if stem1.forms_pseudoknot_with(stem2):
                pseudoknots.append(Pseudoknot(stem1, stem2))
        return pseudoknots


def generate_test_function(bpseq_path):
    bpseq = BPSEQ.from_file(bpseq_path)
    with open(bpseq_path.replace('.bpseq', '-stems.json')) as infile:
        stems = json.load(infile, cls=Decoder)
    with open(bpseq_path.replace('.bpseq', '-hairpins.json')) as infile:
        hairpins = json.load(infile, cls=Decoder)
    with open(bpseq_path.replace('.bpseq', '-pseudoknots.json')) as infile:
        pseudoknots = json.load(infile, cls=Decoder)

    test_functions = []
    if stems:
        def test_function_stems():
            assert bpseq.stems() == stems

        test_function_stems.__name__ = 'test_stems_{}'.format(os.path.basename(bpseq_path))
        test_functions.append(test_function_stems)

    if hairpins:
        def test_function_hairpins():
            #print(hairpins)
            assert bpseq.hairpins() == hairpins

        test_function_hairpins.__name__ = 'test_hairpins_{}'.format(os.path.basename(bpseq_path))
        test_functions.append(test_function_hairpins)

    if pseudoknots:
        def test_function_pseudoknots():
            assert bpseq.pseudoknots() == pseudoknots

        test_function_pseudoknots.__name__ = 'test_pseudoknots_{}'.format(os.path.basename(bpseq_path))
        test_functions.append(test_function_pseudoknots)

    return test_functions


if __name__ == '__main__':
    suite = unittest.TestSuite()
    for bpseq_path in glob.iglob('data/*.bpseq'):
        for test_function in generate_test_function(bpseq_path):
            suite.addTest(unittest.FunctionTestCase(test_function))
    unittest.TextTestRunner().run(suite)
