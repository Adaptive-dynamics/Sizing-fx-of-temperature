HASDOCKER ?= $(shell which docker)

DOC := $(if $(HASDOCKER), docker run --net host --rm --user=$$(id -u):$$(id -g) -v $$PWD:/work -w /work docker.dragonfly.co.nz/dragonverse-17.04:2017-06-28,)
all: poster.pdf

#We need to fall back to manually calling the programs if using old latexmk
poster.pdf:poster.tex
	$(DOC) bash -c "(TEXINPUTS=.///: xelatex $<) && (TEXINPUTS=.///: biber $(<:%.tex=%)) && (TEXINPUTS=.///: xelatex $<)"

#Sometimes we want to blow away all of the R code output as well as just the latex stuff
superclean: clean
	$(DOC) bash -c "rm -rf *.aux *.bcf *.log *.blg *.xml *.bbl *.out *.bib"
clean:
	$(DOC) bash -c "rm -rf  $(OUTPUT_DIR) $(TARGET)"

