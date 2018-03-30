#!/bin/bash

DEFAULT_PPLACER_BUILD=1.1.alpha19

if [[ -z $1 ]]; then
    echo "usage: $(basename $0) <prefix> [<pplacer-build>]"
    echo
    echo "Installs pplacer binaries to <prefix>/bin"
    echo "and accompanying python scripts to wherever pip2 puts things"
    echo "the default value of pplacer-build is $DEFAULT_PPLACER_BUILD"
    exit 1
fi

PREFIX=$(readlink -f $1)
PPLACER_BUILD=${2-$DEFAULT_PPLACER_BUILD}
PPLACER_DIR=pplacer-Linux-v${PPLACER_BUILD}
PPLACER_ZIP=${PPLACER_DIR}.zip
PPLACER_URL=https://github.com/matsen/pplacer/releases/download/v${PPLACER_BUILD}/$PPLACER_ZIP

pplacer_is_installed(){
    $PREFIX/bin/pplacer --version 2> /dev/null | grep -q "$PPLACER_BUILD"
}

if pplacer_is_installed; then
    echo -n "pplacer is already installed: "
    $PREFIX/bin/pplacer --version
else
    mkdir -p $PREFIX/bin && \
    mkdir -p src && \
        (cd src && \
        wget --quiet -nc $PPLACER_URL && \
        unzip -o $PPLACER_ZIP && \
        cp $PPLACER_DIR/{pplacer,guppy,rppr} $PREFIX/bin && \
        pip2 install -U $PPLACER_DIR/scripts)
    # confirm that we have installed the requested build
    if ! pplacer_is_installed; then
        echo -n "Error: you requested pplacer build $PPLACER_BUILD "
        echo "but $($PREFIX/bin/pplacer --version) is installed."
        echo "Try removing src/$PPLACER_ZIP first."
        exit 1
    fi
fi
