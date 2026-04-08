#!/usr/bin/env python
import csv
import sys
import unittest
from io import StringIO

from sqlalchemy import create_engine

from taxtastic.subcommands import append
from taxtastic.taxonomy import Taxonomy

from . import config
from .config import TestBase

dbname = config.ncbi_master_db


class Args:
    schema = None
    verbosity = 0
    tax_id_column = 'tax_id'


class TestAppendAction(TestBase):

    def setUp(self):
        self.engine = create_engine('sqlite:///%s' % dbname)
        self.tax = Taxonomy(self.engine)

    def tearDown(self):
        self.engine.dispose()

    def _make_args(self, rows, columns):
        buf = StringIO()
        writer = csv.DictWriter(buf, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
        buf.seek(0)

        args = Args()
        args.url = 'sqlite:///%s' % dbname
        args.infile = buf
        args.columns = columns
        args.outfile = StringIO()
        return args

    def test_appends_genus_and_species_name(self):
        """Normal case: valid tax_ids get genus and species_name appended."""
        rows = [
            {'seqname': 'seq1', 'tax_id': '1280'},  # Staphylococcus aureus
            {'seqname': 'seq2', 'tax_id': '1378'},  # Gemella
        ]
        args = self._make_args(rows, ['genus', 'species_name'])
        append.action(args)

        args.outfile.seek(0)
        result = list(csv.DictReader(args.outfile))

        self.assertEqual(len(result), 2)
        self.assertIn('genus', result[0])
        self.assertIn('species_name', result[0])

        # Staphylococcus aureus (1280) should have a genus tax_id
        self.assertTrue(result[0]['genus'])
        self.assertEqual(
            result[0]['species_name'], 'Staphylococcus aureus')

    def test_existing_columns_not_duplicated(self):
        """Columns already in the CSV are not added twice to fieldnames."""
        rows = [{'tax_id': '1280', 'genus': 'existing_value'}]
        args = self._make_args(rows, ['genus'])
        append.action(args)

        args.outfile.seek(0)
        reader = csv.DictReader(args.outfile)
        self.assertEqual(reader.fieldnames.count('genus'), 1)

    def test_missing_tax_id_produces_empty_columns(self):
        """Tax_ids not found in the database produce empty column values."""
        rows = [{'tax_id': 'invalid_tax_id_999'}]
        args = self._make_args(rows, ['genus', 'species_name'])
        append.action(args)

        args.outfile.seek(0)
        result = list(csv.DictReader(args.outfile))

        self.assertEqual(len(result), 1)
        self.assertEqual(result[0]['genus'], '')
        self.assertEqual(result[0]['species_name'], '')

    def test_output_preserves_input_columns(self):
        """All original columns are preserved in the output."""
        rows = [
            {'seqname': 'seq1', 'description': 'a desc', 'tax_id': '1280'},
        ]
        args = self._make_args(rows, ['species_name'])
        append.action(args)

        args.outfile.seek(0)
        result = list(csv.DictReader(args.outfile))

        self.assertEqual(result[0]['seqname'], 'seq1')
        self.assertEqual(result[0]['description'], 'a desc')
        self.assertEqual(result[0]['tax_id'], '1280')
