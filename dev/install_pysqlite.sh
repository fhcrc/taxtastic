#!/bin/bash

set -e

sqlite_min_version=3.8.3
pysqlite_version=2.8.3

sqlite_ok=$(python <<EOF
from distutils.version import LooseVersion as v
try:
    from pysqlite2 import dbapi2 as sqlite3
    module='pysqlite2'
except ImportError:
    import sqlite3
    module='sqlite3'

if v(sqlite3.sqlite_version) > v("$sqlite_min_version"):
    print '{}-{}'.format(module, sqlite3.sqlite_version)

EOF
)

echo $sqlite_ok

if [[ $sqlite_ok ]]; then
    echo "sqlite library version is > $sqlite_min_version - exiting"
    exit 0
fi

rm -rf pysqlite-$pysqlite_version
pip2 install -U pip
pip2 download pysqlite==$pysqlite_version
tar -xf pysqlite-$pysqlite_version.tar.gz
cd pysqlite-$pysqlite_version

# see https://github.com/ghaering/pysqlite/issues/95
cat > setup_patch <<EOF
146a147,149
>             ext.include_dirs.append(".")
>             build_ext.build_extension(self, ext)
>             return
EOF

patch setup.py setup_patch

wget --quiet https://sqlite.org/2017/sqlite-amalgamation-3200100.zip
unzip sqlite-amalgamation-3200100.zip
mv sqlite-amalgamation-3200100/* .

rm -rf build && python setup.py build_static && python setup.py install

python -c 'from pysqlite2 import dbapi2; print dbapi2.sqlite_version'
