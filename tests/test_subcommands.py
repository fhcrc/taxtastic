import contextlib
import unittest
import tempfile
import shutil
import copy
import os
import os.path
import csv
import sys

import sqlalchemy as sa

from taxtastic import refpkg
from taxtastic.subcommands import (
    update, create, strip, rollback, rollforward,
    taxtable, check, add_to_taxtable)
from taxtastic.scripts.taxit import main
from taxtastic.taxonomy import Taxonomy


from . import config
from .config import OutputRedirectMixin, TestBase, data_path


class TestUpdate(OutputRedirectMixin, unittest.TestCase):

    def setUp(self):
        super(TestUpdate, self).setUp()

        class _Args(object):
            refpkg = None
            changes = []
            metadata = False
            stats_type = None
            frequency_type = None
        self.args = _Args()

    def test_action(self):
        with config.tempdir() as scratch:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path, create=True)
            test_file = data_path('bv_refdata.csv')

            self.args.refpkg = pkg_path
            self.args.changes = ['meep=' + test_file, 'hilda=' + test_file]

            update.action(self.args)
            r._sync_from_disk()
            self.assertEqual(r.contents['files']['meep'], 'bv_refdata.csv')

            # Second file should have been assigned a non-clashing name
            h = r.contents['files']['hilda']
            self.assertNotEqual(h, 'bv_refdata.csv')
            self.assertTrue(h.startswith('bv_refdata'))
            self.assertTrue(h.endswith('.csv'))

            self.assertTrue(os.path.exists(r.resource_path('hilda')))

    def test_metadata_action(self):
        with config.tempdir() as scratch:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path, create=True)
            self.args.changes = ['meep=boris', 'hilda=vrrp']
            self.args.metadata = True
            self.args.refpkg = pkg_path
            update.action(self.args)
            r._sync_from_disk()
            self.assertEqual(r.metadata('meep'), 'boris')
            self.assertEqual(r.metadata('hilda'), 'vrrp')

    def test_update_stats_action(self):
        with config.tempdir() as scratch:
            pkg_path = os.path.join(scratch, 'test.refpkg')
            r = refpkg.Refpkg(pkg_path, create=True)

            args = self.args
            stats_path = os.path.join(config.datadir, 'phyml_aa_stats.txt')

            args.refpkg = pkg_path
            args.changes = ['tree_stats=' + stats_path]
            args.frequency_type = 'empirical'

            update.action(args)

            r._sync_from_disk()

            self.assertIn('tree_stats', r.contents['files'])
            self.assertIn('phylo_model', r.contents['files'])
            self.assertTrue(r.contents['files'][
                            'phylo_model'].endswith('.json'))


class TestCreate(OutputRedirectMixin, unittest.TestCase):

    def setUp(self):
        super(TestCreate, self).setUp()

        class _Args(object):
            clobber = True
            locus = 'Nowhere'
            description = 'A description'
            author = 'Boris the Mad Baboon'
            package_version = '0.3'
            tree_stats = None
            stats_type = None
            aln_fasta = None
            aln_sto = None
            phylo_model = None
            seq_info = None
            mask = None
            model = None

            profile = None
            readme = None
            tree = None
            taxonomy = None
            reroot = False
            rppr = 'rppr'
            frequency_type = None

            def __init__(self, scratch):
                self.package_name = os.path.join(scratch, 'test.refpkg')
        self._Args = _Args

    def test_create(self):
        with config.tempdir() as scratch:
            args = self._Args(scratch)
            create.action(args)
            r = refpkg.Refpkg(args.package_name, create=False)
            self.assertEqual(r.metadata('locus'), 'Nowhere')
            self.assertEqual(r.metadata('description'), 'A description')
            self.assertEqual(r.metadata('author'), 'Boris the Mad Baboon')
            self.assertEqual(r.metadata('package_version'), '0.3')
            self.assertEqual(r.metadata('format_version'), '1.1')
            self.assertEqual(r.contents['rollback'], None)
            args2 = self._Args(scratch)
            args2.package_name = os.path.join(scratch, 'test.refpkg')
            args2.clobber = True
            self.assertEqual(0, create.action(args2))

    def _test_create_phylo_model(self, stats_path, stats_type=None,
                                 frequency_type=None):
        with config.tempdir() as scratch:
            args = self._Args(scratch)
            args.tree_stats = stats_path
            args.stats_type = stats_type
            args.frequency_type = frequency_type
            create.action(args)

            r = refpkg.Refpkg(args.package_name, create=False)
            self.assertIn('phylo_model', r.contents['files'])

    def test_create_phyml_aa(self):
        stats_path = os.path.join(config.datadir, 'phyml_aa_stats.txt')
        self._test_create_phylo_model(stats_path, frequency_type='model')
        self._test_create_phylo_model(
            stats_path, 'PhyML', frequency_type='model')
        self._test_create_phylo_model(
            stats_path, 'PhyML', frequency_type='empirical')
        self.assertRaises(ValueError, self._test_create_phylo_model,
                          stats_path, 'FastTree', frequency_type='empirical')
        self.assertRaises(ValueError, self._test_create_phylo_model,
                          stats_path, 'garli', frequency_type='empirical')


class TestStrip(OutputRedirectMixin, unittest.TestCase):

    def test_strip(self):
        with config.tempdir() as scratch:
            rpkg = os.path.join(scratch, 'tostrip.refpkg')
            shutil.copytree(data_path(
                'lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg, create=False)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')

            class _Args(object):
                refpkg = rpkg
            strip.action(_Args())

            r._sync_from_disk()
            self.assertEqual(r.contents['rollback'], None)
            self.assertEqual(r.contents['rollforward'], None)


class TestRollback(OutputRedirectMixin, unittest.TestCase):
    maxDiff = None

    def test_rollback(self):
        with config.tempdir() as scratch:
            rpkg = os.path.join(scratch, 'tostrip.refpkg')
            shutil.copytree(data_path(
                'lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg, create=False)
            original_contents = copy.deepcopy(r.contents)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')
            updated_contents = copy.deepcopy(r.contents)

            class _Args(object):
                refpkg = rpkg

                def __init__(self, n):
                    self.n = n

            self.assertEqual(rollback.action(_Args(3)), 1)
            r._sync_from_disk()
            self.assertEqual(r.contents, updated_contents)

            self.assertEqual(rollback.action(_Args(2)), 0)
            r._sync_from_disk()
            self.assertEqual(r.contents['metadata'],
                             original_contents['metadata'])
            self.assertEqual(r.contents['rollback'], None)
            self.assertNotEqual(r.contents['rollforward'], None)


class TestRollforward(TestBase):
    maxDiff = None

    def setUp(self):
        self.suppress_stderr()
        self.suppress_stdout()

    def test_rollforward(self):
        with config.tempdir() as scratch:
            rpkg = os.path.join(scratch, 'tostrip.refpkg')
            shutil.copytree(data_path(
                'lactobacillus2-0.2.refpkg'), rpkg)
            r = refpkg.Refpkg(rpkg, create=False)
            original_contents = copy.deepcopy(r.contents)
            r.update_metadata('boris', 'hilda')
            r.update_metadata('meep', 'natasha')
            updated_contents = copy.deepcopy(r.contents)
            r.rollback()
            r.rollback()

            class _Args(object):
                refpkg = rpkg

                def __init__(self, n):
                    self.n = n

            self.assertEqual(rollforward.action(_Args(3)), 1)
            r._sync_from_disk()
            self.assertEqual(r.contents['metadata'],
                             original_contents['metadata'])

            self.assertEqual(rollforward.action(_Args(2)), 0)
            r._sync_from_disk()
            self.assertEqual(r.contents['metadata'],
                             updated_contents['metadata'])
            self.assertEqual(r.contents['rollforward'], None)
            self.assertNotEqual(r.contents['rollback'], None)


@contextlib.contextmanager
def scratch_file(unlink=True):
    """Create a temporary file and return its name.

    At the start of the with block a secure, temporary file is created
    and its name returned.  At the end of the with block it is
    deleted.
    """
    try:
        (tmp_fd, tmp_name) = tempfile.mkstemp(text=True)
        os.close(tmp_fd)
        yield tmp_name
    except ValueError as v:
        raise v
    else:
        if unlink:
            os.unlink(tmp_name)


class TestTaxtable(TestBase):

    def test_invalid_taxid(self):
        outdir = self.mkoutdir()
        args = ['taxtable', config.ncbi_master_db,
                '--tax-ids', 'horace', '1280',
                '-o', os.path.join(outdir, 'taxonomy.csv')]
        self.assertRaises(ValueError, main, args)

    def test_seqinfo(self):
        outdir = self.mkoutdir()
        args = ['taxtable', config.ncbi_master_db,
                '--seq-info', data_path('simple_seqinfo.csv'),
                '-o', os.path.join(outdir, 'taxonomy.csv')]
        self.assertIsNone(main(args))


class TestAddToTaxtable(TestBase):
    maxDiff = None

    def test_seqinfo(self):
        args = ['add_to_taxtable',
                data_path('minimal_taxonomy.csv'),
                data_path('minimal_add_taxonomy.csv'),
                '-o', os.path.join(self.mkoutdir(), 'taxonomy.csv')]
        self.assertIsNone(main(args))


class TestCheck(OutputRedirectMixin, unittest.TestCase):

    def test_runs(self):
        class _Args(object):
            refpkg = data_path('lactobacillus2-0.2.refpkg')
        self.assertEqual(check.action(_Args()), 0)


class TestUpdateTaxids(TestBase):

    def setUp(self):
        if '-v' not in sys.argv:
            self.suppress_stdout()
            self.suppress_stderr()

        self.outdir = self.mkoutdir()
        self.infile = os.path.join(self.outdir, 'infile.csv')
        self.outfile = os.path.join(self.outdir, 'outfile.csv')
        self.unknowns = os.path.join(self.outdir, 'unknowns.csv')
        self.db = data_path('small_taxonomy.db')

        self.input = [
            ('tax_id', 'tax_name', 'comment'),
            ('1280', '', 'ok'),
            ('1291', 'Staphylococcus staphylolyticus', 'merged with 1287'),
            ('', 'who knows?', 'blank'),
            ('foo', 'unknown', 'completely unknown'),
        ]

        with open(self.infile, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(self.input)

    def get_rows(self, fname):
        with open(fname) as f:
            reader = csv.reader(f)
            return [tuple(row) for row in reader]

    def test01(self):
        self.assertRaises(SystemExit, main, ['update_taxids', '-h'])

    def test02(self):
        # test the test harness itself
        self.assertEqual(self.input, self.get_rows(self.infile))

    def test03(self):
        args = ['update_taxids', self.infile, self.db]
        self.assertRaises(SystemExit, main, args)

    def test04(self):
        args = ['update_taxids', self.infile, self.db,
                '-o', self.outfile,
                '--unknown-action', 'ignore']
        main(args)

        expected = [
            ('tax_id', 'tax_name', 'comment'),
            ('1280', '', 'ok'),
            ('1287', 'Staphylococcus staphylolyticus', 'merged with 1287'),
            ('', 'who knows?', 'blank'),
            ('foo', 'unknown', 'completely unknown'),
        ]

        self.assertEqual(expected, self.get_rows(self.outfile))

    def test05(self):
        args = ['update_taxids', self.infile, self.db,
                '-o', self.outfile,
                '--unknown-action', 'drop']
        main(args)

        expected = [
            ('tax_id', 'tax_name', 'comment'),
            ('1280', '', 'ok'),
            ('1287', 'Staphylococcus staphylolyticus', 'merged with 1287'),
        ]

        self.assertEqual(expected, self.get_rows(self.outfile))

    def test06(self):
        args = ['update_taxids', self.infile, self.db,
                '-o', self.outfile,
                '--unknown-action', 'ignore',
                '--unknowns', self.unknowns]
        main(args)

        unknowns = [
            ('tax_id', 'tax_name', 'comment'),
            ('', 'who knows?', 'blank'),
            ('foo', 'unknown', 'completely unknown'),
        ]

        self.assertEqual(unknowns, self.get_rows(self.unknowns))


class TestAddNode(TestBase):

    def setUp(self):
        if '-v' not in sys.argv:
            self.suppress_stdout()
            self.suppress_stderr()

        self.outdir = self.mkoutdir()
        self.dbname = os.path.join(self.outdir, 'taxonomy.db')
        shutil.copyfile(data_path('small_taxonomy.db'), self.dbname)

    def assertZeroExitStatus(self, val):
        self.assertFalse(bool(val))

    def assertNonZeroExitStatus(self, val):
        self.assertTrue(bool(val))

    def test_new_nodes01(self):
        args = ['add_nodes', self.dbname, data_path('new_nodes_ok.yml')]
        self.assertZeroExitStatus(main(args))

    def test_new_nodes02(self):
        # fails without --source-name
        args = ['add_nodes', self.dbname,
                data_path('new_nodes_ok_nosource.yml')]
        self.assertNonZeroExitStatus(main(args))

    def test_new_nodes03(self):
        args = ['add_nodes', self.dbname,
                data_path('new_nodes_ok_nosource.yml'),
                '--source-name', 'some_source']
        self.assertZeroExitStatus(main(args))

    def test_new_nodes04(self):
        args = ['add_nodes', self.dbname,
                data_path('staph_species_group.yml'),
                '--source-name', 'foo']
        self.assertZeroExitStatus(main(args))

        tax = Taxonomy(sa.create_engine('sqlite:///' + self.dbname))
        with tax.engine.connect() as con:
            result = tax.fetchall(
                sa.text('select source_id from nodes '
                        'where parent_id = :parent_id'),
                parent_id='stapha_sg'
            )

        self.assertEqual(len(result), 5)
        self.assertEqual([row[0] for row in result], [2] * len(result))

    def test_new_nodes05(self):
        args = ['add_nodes', self.dbname, data_path('staph_species_group2.yml')]
        self.assertZeroExitStatus(main(args))

        tax = Taxonomy(sa.create_engine('sqlite:///' + self.dbname))
        with tax.engine.connect() as con:
            result = tax.fetchall(
                sa.text('select source_id from nodes '
                        'where parent_id = :parent_id'),
                parent_id='stapha_sg'
            )

        self.assertEqual(len(result), 5)
        self.assertEqual([row[0] for row in result], [2] * len(result))


class TestExtractNodes(TestBase):

    def setUp(self):
        if '-v' not in sys.argv:
            self.suppress_stdout()
            self.suppress_stderr()

        self.outdir = self.mkoutdir()
        self.dbname = os.path.join(self.outdir, 'taxonomy.db')
        self.outfile = os.path.join(self.outdir, 'extracted.yml')
        shutil.copyfile(data_path('small_taxonomy.db'), self.dbname)

    def test_new_nodes01(self):
        source_name = 'some_source'

        # add some nodes and names
        args = ['add_nodes', self.dbname, data_path('new_nodes_ok_nosource.yml'),
                '--source-name', source_name]
        main(args)
        args = ['extract_nodes', self.dbname, source_name,
                '-o', self.outfile]
        main(args)


class TestLineageTable(TestBase):
    def setUp(self):
        self.outdir = self.mkoutdir()
        self.taxtable = os.path.join(self.outdir, 'taxtable.csv')
        self.seq_info = os.path.join(self.outdir, 'seq_info.csv')
        self.info = [
            ('seqname', 'tax_id', 'species', 'mothur'),
            ('s1', '1280', 'Staphylococcus aureus', 's__Staphylococcus_aureus;'),
            ('s2', '246432', 'Staphylococcus equorum', 's__Staphylococcus_equorum;'),
            ('s3', '29383', 'Staphylococcus equorum', 's__Staphylococcus_equorum;'),
            ('s4', '1279', '', 'g__Staphylococcus;')
        ]

        with open(self.seq_info, 'w') as f:
            csv.writer(f).writerows(self.info)

        main(['taxtable', config.ncbi_master_db, '-i',
              self.seq_info, '-o', self.taxtable])

    def test_csv_table(self):
        outfile = os.path.join(self.outdir, 'lineages.csv')
        args = ['lineage_table', self.taxtable, self.seq_info,
                '--csv-table', outfile]
        main(args)

        with open(outfile) as f:
            output = list(csv.DictReader(f))
            self.assertEqual(
                [row['seqname'] for row in output],
                [row[0] for row in self.info[1:]]
            )

            self.assertEqual(
                [row['species'] for row in output],
                [row[-2] for row in self.info[1:]]
            )

    def test_lineage_table(self):
        outfile = os.path.join(self.outdir, 'taxonomy.txt')
        args = ['lineage_table', self.taxtable, self.seq_info,
                '--taxonomy-table', outfile]
        main(args)

        with open(outfile) as f:
            output = list(csv.reader(f, delimiter='\t'))

            self.assertEqual(
                [row[0] for row in output],
                [row[0] for row in self.info[1:]]
            )

            for expected, actual in zip(self.info[1:], output):
                self.assertTrue(actual[1].endswith(actual[-1]))
