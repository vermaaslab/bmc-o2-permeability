# Convenience wrappers; every target just calls scripts/make_all.py.
PY ?= python

.PHONY: all quick checks clean
all:      ## regenerate Tables 1-2 and Figures 4-5 from data/ (~10 min on 32 cores)
	$(PY) scripts/make_all.py
quick:    ## redo only the cheap steps; reuse the committed ODE-fit outputs
	$(PY) scripts/make_all.py --skip-fits
checks:   ## everything above plus the slow robustness checks (~30 min)
	$(PY) scripts/make_all.py --checks
clean:    ## remove the git-ignored cache (committed results are untouched)
	rm -rf results/cache scripts/__pycache__
