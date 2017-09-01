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

    print config.data_path('small_taxonomy.db')

    def test01(self):
        """
        Exit on error
        """
        args = [self.seq_info, self.small_taxonomy_db]
        log.info(self.log_info + ' '.join(map(str, args)))
        self.assertRaises(SystemExit, self.main, args)

    def test02(self):
        """
        ignore unknowns
        """

        this_test = sys._getframe().f_code.co_name
        thisdata_path = self.thisdata_path
        ref = os.path.join(thisdata_path, this_test, 'update.csv')
        outdir = self.mkoutdir()
        out = os.path.join(outdir, 'update.csv')

        args = ['--unknown-action', 'ignore',
                '--outfile', out,
                self.seq_info, self.small_taxonomy_db]
        log.info(self.log_info + ' '.join(map(str, args)))
        self.main(args)

        print 'diff -w {} {}'.format(ref, out)
        self.assertTrue(filecmp.cmp(ref, out))

    def test03(self):
        """
        write unknowns to a file
        """

        this_test = sys._getframe().f_code.co_name
        thisdata_path = self.thisdata_path
        ref_info = os.path.join(thisdata_path, this_test, 'update.csv')
        ref_unknowns = os.path.join(thisdata_path, this_test, 'unknowns.csv')
        outdir = self.mkoutdir()
        out_info = os.path.join(outdir, 'update.csv')
        out_unknowns = os.path.join(outdir, 'unknowns.csv')
        args = [
            '--unknown-action', 'ignore',
            '--unknowns', out_unknowns,
            '--outfile', out_info,
            self.seq_info,
            self.small_taxonomy_db]
        log.info(self.log_info + ' '.join(map(str, args)))
        self.main(args)

        print 'diff -w {} {}'.format(ref_unknowns, out_unknowns)
        self.assertTrue(filecmp.cmp(ref_unknowns, out_unknowns))

        print 'diff -w {} {}'.format(ref_info, out_info)
        self.assertTrue(filecmp.cmp(ref_info, out_info))

    # did not implement feature to check names for now
    # def test04(self):
    #     """
    #     ignore unknowns, search with tax name
    #     """

    #     this_test = sys._getframe().f_code.co_name
    #     thisdata_path = self.thisdata_path
    #     ref_info = os.path.join(thisdata_path, this_test, 'update.csv')
    #     outdir = self.mkoutdir()
    #     out_info = os.path.join(outdir, 'update.csv')
    #     args = [
    #         '--unknown-action', 'ignore',
    #         '--name-column', 'tax_name',
    #         '--outfile', out_info,
    #         self.seq_info,
    #         self.small_taxonomy_db]
    #     log.info(self.log_info + ' '.join(map(str, args)))
    #     self.main(args)
    #     self.assertTrue(filecmp.cmp(ref_info, out_info))
