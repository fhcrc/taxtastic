import unittest

from taxtastic import greengenes

class ParseClassesTestCase(unittest.TestCase):

    def _test_parse(self, lineage_str, otu_number, expected):
        actual = greengenes._parse_classes(lineage_str, otu_number)
        self.assertEqual(len(expected), len(actual))
        for e, a in zip(expected, actual):
            self.assertEqual(len(e), len(a))
            e_cat, e_name = e
            a_cat, a_name = a
            self.assertEqual(e_cat, a_cat)
            self.assertEqual(e_name, a_name)

    def test_fully_defined(self):
        s = """k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacteriales;f__Enterobacteriaceae;g__Escherichia;s__Citrobacterkoseri"""
        otu = 471423
        expected = [['root', 'root'], ['kingdom', 'Bacteria'],
                    ['phylum', 'Proteobacteria'],
                    ['class', 'Gammaproteobacteria'],
                    ['order', 'Enterobacteriales'],
                    ['family', 'Enterobacteriaceae'],
                    ['genus', 'Escherichia'],
                    ['species', 'Citrobacterkoseri'],
                    ('otu', '471423')]
        self._test_parse(s, otu, expected)

    def test_genus_in_species(self):
        s = """k__Bacteria;g__Escherichia;s__Escherichiakoseri"""
        otu = 471323
        expected = [['root', 'root'], ['kingdom', 'Bacteria'],
                    ['genus', 'Escherichia'],
                    ['species', 'Escherichia koseri'],
                    ('otu', '471323')]
        self._test_parse(s, otu, expected)
