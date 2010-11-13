RUNIT_DIR=$(shell pwd)

ifdef TARGET
	TARG=-t ${TARGET}
endif

ifdef BIN
	BIN=${BIN}/
endif

all: test

install:
	cd ..; \
	python setup.py install && \
	cd ${RUNIT_DIR}

test:
	${BIN}runtests.py ${TARG} ${FLAGS}

%.py: FORCE
	${BIN}runtests.py $@ ${TARG} ${FLAGS}

FORCE:

clean:
	rm -f *.pyc profile ../test_output/* 

