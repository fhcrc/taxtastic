#!/bin/bash

# run from the top level directory to build the sphinx docs, commit
# and push changes, eg:
#   cd bioscons
#   dev/ghpages-publish.sh

(cd docs && make html) && \
(cd html && git add . && git commit -m "Automated sphinx doc build" && git push) && \
git add html && \
echo "changes to html staged - commit and push manually, eg"
echo "git commit -m \"updated sphinx docs\" && git push [--force]"

