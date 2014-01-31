autopep8:
	autopep8 --aggressive --in-place -r taxtastic tests

test:
	./test -v

clean: clean-pyc
	python setup.py clean

clean-pyc:
	find taxtastic tests -name \*.pyc | xargs rm

.PHONY: test autopep8 clean-pyc clean
