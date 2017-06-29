
SHELL := /bin/bash

HASDOCKER ?= $(shell which docker)

DOC := $(if $(HASDOCKER), docker run --net host --rm -v $$PWD:/work -w /work docker.dragonfly.co.nz/dragonverse-17.04:2017-06-28,)

OUTPUT_DIR=build
TARGET=draft.pdf

KNITR = $(addsuffix .tex, $(basename $(shell find . -iname "*.Rnw")))
KNITR_BASE = $(addprefix $(OUTPUT_DIR)/, $(addsuffix .tex, $(basename $(wildcard *.Rnw))))
KNIT_COMMAND = library(knitr);opts_chunk\$$set(warning=F, message = FALSE,echo=T,results='asis',error=FALSE,fig.lp='fig:',fig.path='images/');

####################################################
#### BUILD PROCESS #################################
####################################################

.SECONDEXPANSION:
# do not remove intermediates
.SECONDARY:

all: $(TARGET) 

$(TARGET): $(OUTPUT_DIR)/$(TARGET)
	cp $(OUTPUT_DIR)/$(TARGET) $(TARGET)

$(OUTPUT_DIR)/$(TARGET): $(TARGET:%.pdf=%.tex) $(OUTPUT_DIR) $(KNITR_BASE) #$(KNITR_BASE)
	$(DOC) bash -c "(TEXINPUTS=.///: xelatex -output-directory=$(OUTPUT_DIR) $<) && (TEXINPUTS=.///: biber --output_directory=$(OUTPUT_DIR) $(<:%.tex=%)) && (TEXINPUTS=.///: xelatex -output-directory=$(OUTPUT_DIR) $<) && (TEXINPUTS=.///: xelatex -output-directory=$(OUTPUT_DIR) $<)"
$(OUTPUT_DIR):
	$(DOC) mkdir -p $(OUTPUT_DIR)

$(OUTPUT_DIR)/%.tex: %.Rnw $(KNITR)
	$(DOC) Rscript --vanilla -e "$(KNIT_COMMAND) knit('$(<F)',output='$(OUTPUT_DIR)/$(@F)')"

%.tex: %.Rnw
	$(DOC) Rscript --vanilla -e "$(KNIT_COMMAND) knit('$(<F)',output='$(@F)')"


#knitr/%_tab.tex: knitr/%_tab.rnw 
#	$(DOC) Rscript --vanilla -e "$(KNIT_COMMAND) knit('$(@D)/$(<F)', output = '$(@D)/$(@F)')"

#knitr/%.Rdata: knitr/set_%.R
#	cd $(@D); Rscript $(<F)
	
#knitr/%.tex: knitr/%.rnw 
#	$(DOC) Rscript --vanilla -e "$(KNIT_COMMAND) knit('$(@D)/$(<F)', output = '$(@D)/$(@F)')"

superclean: clean
	$(DOC) find knitr -iname "*.tex" -delete
	$(DOC) find knitr -iname "*.pdf" -delete
	$(DOC) find knitr -iname "*.png" -delete
	$(DOC) find knitr -iname "images" -delete

clean:
	$(DOC) rm -rf  $(OUTPUT_DIR) $(TARGET) 


