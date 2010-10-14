import ConfigParser
import logging
import os
import time
import shutil
import hashlib

log = logging

manifest_name = 'CONTENTS.txt'

package_contents = {
    'metadata':['create_date','author','description'],
    'files':['tree_file','tree_stats','aln_fasta','aln_sto',
             'profile', 'seq_info','taxonomy','mask'],
    'md5':[]
    }

def write_config(fname, optdict, sections):
    """
    * fname - name of config file
    * optdict - a dict with keys d[(section, opt)] = val
    * sections - a list of tuples, eg
      (('section1',['opt1','opt2',...]),('section2',['opt3','opt4',...]))
      where 'opt*' are attributes of 'options'
    """

    config = ConfigParser.SafeConfigParser()

    for section_name, opts in sections.items():
        config.add_section(section_name)
        for opt in opts:
            val = optdict.get((section_name,opt),'')
            config.set(section_name, opt, str(val))

    cfile = open(fname,'w+')
    config.write(cfile)
    log.info('writing %s' % fname)

    cfile.seek(0)
    # log.info('\n' + cfile.read())
    print cfile.read()
    cfile.close()

def create(pkg_dir, options, manifest_name=manifest_name, package_contents=package_contents):

    os.mkdir(pkg_dir)
    manifest = os.path.join(pkg_dir, manifest_name)

    optdict = {}
    optdict[('metadata','create_date')] = time.strftime('%Y-%m-%d %H:%M:%S')

    # copy files into the package directory
    for fname in package_contents['files']:
        pth = getattr(options, fname)

        if pth:
            shutil.copy(pth, pkg_dir)
            optdict[('files',fname)] = os.path.split(pth)[1]
            optdict[('md5',fname)] = hashlib.md5(open(pth).read()).hexdigest()
            package_contents['md5'].append(fname)

    write_config(fname=manifest, optdict=optdict, sections=package_contents)


