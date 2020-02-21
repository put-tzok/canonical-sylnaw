#! /usr/bin/python
import os.path
import sys
import unittest

PDB_IDS = ('1A9N', '1BGZ', '1DRZ', '1EUY', '1FOQ', '1HQ1', '1HVU', '1IVS', '1L2X', '1LNG', '1MV6', '2A43', '2EES',
           '2TPK', '3FU2', '3K1V')


class BP:
    def __init__(self, path=None):
        self.__pairs = []
        if path:
            with open(path) as bpfile:
                for line in bpfile:
                    line = line.strip().split()
                    self.add_pair(int(line[0]), int(line[1]))

    def convert_from_fr3d(self, fr3d_path: str):
        # TODO: implement
        pass

    def add_pair(self, resi, resj):
        self.__pairs.append((resi, resj))

    def __eq__(self, other):
        return frozenset(self.__pairs) == frozenset(other.__pairs)

    def __sub__(self, other):
        return sorted(frozenset(self.__pairs) - frozenset(other.__pairs))


def generate_test_function(pdb_id):
    def test_function():
        path_expected = os.path.join('data', '{}.bp'.format(pdb_id))
        bp_expected = BP(path_expected)

        path_input = os.path.join('data', '{}.fr3d'.format(pdb_id))
        bp_input = BP()
        bp_input.convert_from_fr3d(path_input)

        assert bp_expected == bp_input, '\nProblem found for: {}\nMissing = {}\nRedundant = {}\n'.format(pdb_id,
                                                                                                         bp_expected - bp_input,
                                                                                                         bp_input - bp_expected)

    test_function.__name__ = 'test_{}'.format(pdb_id)
    return test_function


if __name__ == '__main__':
    suite = unittest.TestSuite()
    for pdb_id in PDB_IDS:
        suite.addTest(unittest.FunctionTestCase(generate_test_function(pdb_id)))
    unittest.TextTestRunner().run(suite)
