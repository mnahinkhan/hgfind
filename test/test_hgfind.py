#!/usr/bin/env python3
"""
Tests the hgfind module for correctness.

"""
import unittest

from src.hgfind.hgfind import WrongGeneName, hgfind


class TestGeneToCoord(unittest.TestCase):
    """
    Tests to ensure that hgfind functions correctly

    """

    def test_normal_success(self):
        """
        Checks whether normal gene names are processed correctly

        """

        result = hgfind("HNRNPC")
        self.assertEqual(result["official_name"], "HNRNPC")
        self.assertEqual(result["chr_n"], 14)
        self.assertEqual(result["start_coord"], 21209136)
        self.assertEqual(result["end_coord"], 21269494)

    def test_synonym_success(self):
        """
        Checks that gene synonyms are processed correctly

        """
        result = hgfind("AUF1")
        self.assertEqual(result["official_name"], "HNRNPD")
        self.assertEqual(result["chr_n"], 4)
        self.assertEqual(result["start_coord"], 82352498)
        self.assertEqual(result["end_coord"], 82374503)

    def test_case_sensitivity(self):
        """
        Checks that gene names are not case sensitive

        """

        self.assertEqual(hgfind("auf1"), hgfind("AUF1"))

    def test_gibberish_fail(self):
        """
        Checks that gibberish fails as expected

        """

        self.assertRaises(WrongGeneName, hgfind, "gsjfg")
        self.assertRaises(WrongGeneName, hgfind, "4:45-243")


if __name__ == "__main__":
    unittest.main()
