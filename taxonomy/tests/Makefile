RUNIT_DIR=$(shell pwd)

ifdef TARGET
	TARG=-t ${TARGET}
endif

all: test

install:
	cd ..; \
	python setup.py install && \
	cd ${RUNIT_DIR}

test:
	runtests.py ${TARG} ${FLAGS}

%.py: FORCE
	runtests.py $@ ${TARG} ${FLAGS}

FORCE:

clean:
	rm -f *.pyc profile ../test_output/* 

