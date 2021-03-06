SHELL = /bin/bash -e 

all: build install

build:
	python setup.py build --executable="/usr/bin/env python"

bdist:
	python setup.py build --executable="/usr/bin/env python"
	python setup.py bdist --formats=egg

install:
	python setup.py install --prefix="/galaxy/home/ksahlin/local_python_modules"

develop:
	python setup.py develop

test:
	# Unit tests
	find tests/unit -name "*.py" | xargs nosetests
	# End-to-end tests
	find tests/cram -name "*.t" | xargs cram


doc:
	sphinx-apidoc -T -f -o doc src/ && cd doc && make html
docs: doc

clean: doc-clean
	rm -rf dist/ build/ *.egg-info
	rm -rf doc/_build
	find . -name "*.pyc" | xargs rm -f
	rm -rf dist/
	rm -f nostests.xml
doc-clean:
	rm -f doc/*.html

pip-install:
	@which pip > /dev/null
	@pip freeze|grep 'pbtranscript=='>/dev/null \
      && pip uninstall -y pbtools.pbtranscript \
      || pip uninstall -y pbtranscript \
      || true
	@python setup.py build
	@pip install --no-index \
          --install-option="--install-scripts=$(PREFIX)/bin" \
          ./

.PHONY: all build bdist install develop test doc clean pip-install
