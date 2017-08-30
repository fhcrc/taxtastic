#!/bin/bash -v

wget -nc https://pypi.python.org/packages/42/02/981b6703e3c83c5b25a829c6e77aad059f9481b0bbacb47e6e8ca12bd731/pysqlite-2.8.3.tar.gz

rm -rf pysqlite-2.8.3
tar -xf pysqlite-2.8.3.tar.gz
cd pysqlite-2.8.3

cat > setup_patch <<EOF
146a147,149
>             ext.include_dirs.append(".")
>             build_ext.build_extension(self, ext)
>             return
EOF

patch setup.py setup_patch

wget -nc https://sqlite.org/2017/sqlite-amalgamation-3200100.zip
unzip sqlite-amalgamation-3200100.zip
mv sqlite-amalgamation-3200100/* .

rm -rf build && python setup.py build_static && python setup.py install

python -c 'from pysqlite2 import dbapi2; print dbapi2.sqlite_version'
