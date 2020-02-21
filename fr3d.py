#! /usr/bin/python
import os.path
import sys
import unittest

PDB_IDS = ('1A9N', '1BGZ', '1DRZ', '1EUY', '1FOQ', '1HQ1', '1HVU', '1IVS', '1L2X', '1LNG', '1MV6', '2A43', '2EES',
           '2TPK', '3FU2', '3K1V');


def convert(filecontent):
    # TODO: implement
    pass


def readfile(path):
    with open(path) as fp:
        return fp.read()


def generate_test_function(pdb_id):
    def test_function():
        path_input = os.path.join('data', '{}.fr3d'.format(pdb_id))
        path_expected = os.path.join('data', '{}.bp'.format(pdb_id))
        content_input = readfile(path_input)
        content_expected = readfile(path_expected)
        assert content_expected == convert(content_input)
    test_function.__name__ = 'test_{}'.format(pdb_id)
    return test_function


if __name__ == '__main__':
    suite = unittest.TestSuite()
    for pdb_id in PDB_IDS:
        suite.addTest(unittest.FunctionTestCase(generate_test_function(pdb_id)))
    unittest.TextTestRunner().run(suite)