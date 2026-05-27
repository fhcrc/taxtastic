#!/usr/bin/env python

import re
import os
from os import path
import logging

import sqlalchemy as sa

import taxtastic
import taxtastic.ncbi
from taxtastic.ncbi import read_names, read_archive

from . import config
from .config import TestBase

log = logging

outputdir = config.outputdir
datadir = config.datadir
ncbi_master_db = config.ncbi_master_db
ncbi_data = config.ncbi_data


class TestDbconnect(TestBase):

    def test01(self):
        engine = sa.create_engine('sqlite:///' + ncbi_master_db)
        taxtastic.ncbi.db_connect(engine)
        with engine.begin() as con:
            result = con.execute(sa.text(
                'select name from sqlite_master where type = "table"'))
            tables = set(i[0] for i in result)
            self.assertTrue(
                set(['nodes', 'names', 'merged', 'source']).issubset(tables))


class TestLoadData(TestBase):

    def setUp(self):
        outdir = self.mkoutdir()
        self.db_path = os.path.join(outdir, 'taxonomy.db')
        self.url = 'sqlite:///' + self.db_path

    def test01(self):
        # we should be starting from scratch
        self.assertFalse(path.isfile(self.db_path))
        engine = sa.create_engine(self.url)

        taxtastic.ncbi.db_connect(engine)
        self.assertTrue(path.isfile(self.db_path))


class TestReadNames(TestBase):

    def setUp(self):
        self.zipfile = ncbi_data

    def test02(self):
        """
        is_classified always None
        """

        rows = read_names(rows=read_archive(self.zipfile, 'names.dmp'))
        headers = next(rows)
        is_classified = headers.index('is_classified')
        self.assertEqual(
            set(row[is_classified] for row in rows), set([None]))


class TestUnclassifiedRegex(TestBase):
    """
    Test the heuristic used to determine if a taxonomic name is meaningful.
    """

    unclassified_regex_examples = [
        (r'-like\b', 'WA-like'),
        (r'\bactinomycete\b', 'actinomycete A4'),
        (r'\bcrenarchaeote\b', 'crenarchaeote OlA-6'),
        (r'\bculture\b', 'mixed culture'),
        (r'\bchimeric\b', 'chimeric sequence SJA-7'),
        (r'\bcyanobiont\b', 'Nostocaceae cyanobiont AE1'),
        (r'degrading', '2,4-D degrading bacterium M1'),
        (r'\beuryarchaeote\b', 'euryarchaeote D4.75-4'),
        (r'disease', 'Fiji disease virus'),
        (r'\b[cC]lone', 'soil clone WD1'),
        (r'\bmethanogen(ic)?\b', 'methanogen 5c'),
        (r'\bplanktonic\b', 'planktonic crenarchaeote'),
        (r'\bplanctomycete\b', 'planctomycete A-2'),
        (r'\bsymbiote\b', 'Ornithodoros moubata symbiote B'),
        (r'\btransconjugant\b', '2,4-D degrading transconjugant WD2'),
        (r'^(?!root$)[a-z]', 'liara'),
        # No matching primary scientific name was found in ncbi_taxonomy.db
        # for r'^\W+\s+[a-zA-Z]*\d'.
        (r'\d\d', 'GB10C'),
        (r'atypical', 'Limbodessus atypicalis'),
        (r'\bcf\.', 'Doto cf. divae'),
        (r'acidophile', 'acidophile enrichment culture'),
        (r'\bactinobacterium\b', 'actinobacterium D'),
        (r'aerobic', 'aerobic bacillus'),
        (r'.+\b[Al]g(um|a)\b', 'Prosthecochloris sp. L21-Aga-LIT'),
        (r'\b[Bb]acteri(um|al)\b', 'bacterium'),
        (r'.+\b[Bb]acteria\b', 'cf. Bacteria SAR-1'),
        (r'Barophile', 'Barophile WHB 46'),
        (r'cyanobacterium', 'cyanobacterium G1'),
        (r'Chloroplast', 'Chloroplast expression vector pCEV1'),
        (r'Cloning', 'Cloning vector AC'),
        (r'\bclone\b', 'soil clone WD1'),
        (r'cluster', 'Cyclustera'),
        (r'^diazotroph', 'diazotroph DU1'),
        (r'\bcoccus\b', 'Dactylopius coccus'),
        (r'archaeon', 'archaeon'),
        (r'-containing', 'HSV-tk-containing vector pGTK3'),
        (r'epibiont', 'epibiont metagenome'),
        (r'environmental samples', 'environmental samples <bats>'),
        (r'eubacterium', 'Desulfoeubacterium'),
        (r'halophilic', 'halophilic archaeon'),
        (r'hydrothermal\b', 'hydrothermal vent metagenome'),
        (r'isolate', 'Biserrula isolate'),
        (r'\bmarine\b', 'marine bacterium'),
        (r'methanotroph', 'methanotroph B8'),
        (r'microorganism', 'uncultured microorganism'),
        (r'mollicute', "'Candidatus Bacilliplasma' mollicute"),
        (r'pathogen', 'fungal pathogen UWFP-388'),
        (r'[Pp]hytoplasma', 'Phytoplasma sp.'),
        (r'proteobacterium', 'proteobacterium 1'),
        (r'putative', 'uncultured putative methanotroph'),
        (r'\bsp\.', 'Aix sp.'),
        (r'species', 'Eucara species group'),
        (r'spirochete', 'wall-less spirochete'),
        (r'str\.', 'archaeal str. vp21'),
        (r'strain', 'slope strain DI4'),
        (r'symbiont', 'ant endosymbionts'),
        (r'\b[Tt]axon\b', 'Bisgaard Taxon 2'),
        (r'unicellular', 'Gesasha unicellularis'),
        (r'uncultured', 'uncultured Ulva'),
        (r'unclassified', 'unclassified Aa'),
        (r'unidentified', 'unidentified'),
        (r'unknown', 'HIV-1 unknown group'),
        (r'vector\b', 'YTT vector A'),
        (r'vent\b', 'Ripavent virus'),
    ]

    def setUp(self):
        self.pieces = taxtastic.ncbi.UNCLASSIFIED_REGEX_COMPONENTS
        self.regexes = [re.compile(piece) for piece in self.pieces]
        with open(config.data_path('type_strain_names.txt')) as fp:
            self.type_strain_names = [i.rstrip() for i in fp]

    def test_no_type_strains_match(self):
        for strain_name in self.type_strain_names:
            for regex in self.regexes:
                m = regex.search(strain_name)
                if m:
                    self.fail('"{0}" matches "{1}"'.format(
                        strain_name, regex.pattern))

    def test_cf_prefix_matches(self):
        self.assertRegex(
            'cf. Pedicellasteridae sp. CLM-2010-1',
            taxtastic.ncbi.UNCLASSIFIED_REGEX)

    def test_cf_middle_matches(self):
        self.assertRegex(
            'Plasmodium cf. ovale',
            taxtastic.ncbi.UNCLASSIFIED_REGEX)

    def test_scientific_name_examples_match_components(self):
        for pattern, tax_name in self.unclassified_regex_examples:
            with self.subTest(pattern=pattern, tax_name=tax_name):
                self.assertRegex(tax_name, re.compile(pattern))

# def generate_test_unclassified_regex():
    #"""
    # Generate a test class verifying that none of the type strains in
    # type_strain_names.txt match the unclassified regex.
    #"""
    # def generate_test(strain_name):
        # def do_test(self):
            # for regex in self.regexes:
                #m = regex.search(strain_name)
                # if m:
                    #self.fail('"{0}" matches "{1}"'.format(strain_name, regex.pattern))
        # return do_test

    # class TestUnclassifiedRegex(TestBase):
        # def setUp(self):
            #self.pieces = taxtastic.ncbi.UNCLASSIFIED_REGEX_COMPONENTS
            #self.regexes = [re.compile(piece) for piece in self.pieces]

    # with open(config.data_path('type_strain_names.txt')) as fp:
        #type_strain_names = [i.rstrip() for i in fp]

    # for s in type_strain_names:
        #test_fn = generate_test(s)
        #func_name = 'test_{0}_no_match'.format(re.sub(r'[ -.]', '_', s).lower())
        #test_fn.__name__ = func_name
        #setattr(TestUnclassifiedRegex, func_name, test_fn)

    # return TestUnclassifiedRegex

#TestUnclassifiedRegex = generate_test_unclassified_regex()
#del generate_test_unclassified_regex
