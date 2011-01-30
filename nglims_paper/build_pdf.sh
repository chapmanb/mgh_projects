#!/bin/sh
pdflatex nglims_galaxy.tex
bibtex nglims_galaxy
pdflatex nglims_galaxy.tex
pdflatex nglims_galaxy.tex
