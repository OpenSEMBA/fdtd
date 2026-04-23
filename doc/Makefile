SPHINXBUILD ?= python3 -m sphinx
SOURCEDIR = .
BUILDDIR = _build

.PHONY: help Makefile

help:
	$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)"

%: Makefile
	$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)"
