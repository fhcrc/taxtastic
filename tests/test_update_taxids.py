"""
Test update_taxids
"""

import config
import filecmp
import logging
import os
import sys

from taxtastic.scripts import taxit

log = logging.getLogger(__name__)


class TestUpdateTaxids(config.TestBase):

    def main(self, arguments):
        taxit.main(['update_taxids'] + [str(a) for a in arguments])

    log_info = 'taxit update_taxids '

    thisdata_path = config.data_path('update_taxids', 'TestUpdateTaxids')

    seq_info = config.data_path(thisdata_path, 'seq_info.csv')
    small_taxonomy_db = 'sqlite:///' + config.data_path('small_taxonomy.db')

    def test01(self):
        """
        Minimal inputs
        """
        args = [self.seq_info, self.small_taxonomy_db]
        log.info(self.log_info + ' '.join(map(str, args)))
        # ValueError: Unknown or missing tax_ids present
        self.assertRaises(ValueError, self.main, args)

    def test02(self):
        """
        --ignore-unknowns
        """
        this_test = sys._getframe().f_code.co_name
        thisdata_path = self.thisdata_path
        ref = os.path.join(thisdata_path, this_test, 'update.csv')
        outdir = self.mkoutdir()
        out = os.path.join(outdir, 'update.csv')
        args = ['--ignore-unknowns', '--out', out, self.seq_info, self.small_taxonomy_db]
        log.info(self.log_info + ' '.join(map(str, args)))
        self.main(args)
        self.assertTrue(filecmp.cmp(out, ref))

    def test03(self):
        """
        --unknowns unknowns.csv
        """

        this_test = sys._getframe().f_code.co_name
        thisdata_path = self.thisdata_path
        ref_info = os.path.join(thisdata_path, this_test, 'update.csv')
        ref_unknowns = os.path.join(thisdata_path, this_test, 'unknowns.csv')
        outdir = self.mkoutdir()
        out_info = os.path.join(outdir, 'update.csv')
        out_unknowns = os.path.join(outdir, 'unknowns.csv')
        args = [
            '--unknowns', out_unknowns,
            '--out', out_info,
            self.seq_info,
            self.small_taxonomy_db]
        log.info(self.log_info + ' '.join(map(str, args)))
        self.main(args)
        self.assertTrue(filecmp.cmp(out_info, ref_info))
        self.assertTrue(filecmp.cmp(out_unknowns, ref_unknowns))

    def test04(self):
        """
        --ignore-unknowns --name-columns tax_name
        """

        this_test = sys._getframe().f_code.co_name
        thisdata_path = self.thisdata_path
        ref_info = os.path.join(thisdata_path, this_test, 'update.csv')
        outdir = self.mkoutdir()
        out_info = os.path.join(outdir, 'update.csv')
        args = [
            '--ignore-unknowns',
            '--name-column', 'tax_name',
            '--out', out_info,
            self.seq_info,
            self.small_taxonomy_db]
        log.info(self.log_info + ' '.join(map(str, args)))
        self.main(args)
        self.assertTrue(filecmp.cmp(out_info, ref_info))
