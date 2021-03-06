
SHELL := /bin/bash

HASDOCKER ?= $(shell which docker)

DOC := $(if $(HASDOCKER), docker run --net host --rm --user=$$(id -u):$$(id -g) -v $$PWD:/work -w /work docker.dragonfly.co.nz/dragonverse-18.04,)

OUTPUT_DIR=build
TARGET=draft_model_only.pdf
PRESENTATION=model_meeting.pdf

KNITR = $(addsuffix .tex, $(basename $(shell find ./figures/ -iname "*.Rnw")))
#KNITR_BASE = $(addprefix $(OUTPUT_DIR)/, $(addsuffix .tex, $(basename $(wildcard *.Rnw))))
KNIT_COMMAND = library(knitr);opts_chunk\$$set(warning=F, message = FALSE,echo=T,results='asis',error=FALSE,fig.lp='fig:',fig.path='images/');

####################################################
#### BUILD PROCESS #################################
####################################################

.SECONDEXPANSION:
# do not remove intermediates
.SECONDARY:

all: $(TARGET)

$(TARGET): $(OUTPUT_DIR)/$(TARGET) $(OUTPUT_DIR) $(KNITR)
	cp $(OUTPUT_DIR)/$(TARGET) $(TARGET)

$(PRESENTATION): $(OUTPUT_DIR)/$(PRESENTATION)
	$(DOC) bash -c "cp $(OUTPUT_DIR)/$(PRESENTATION) $(PRESENTATION)"

$(OUTPUT_DIR)/$(PRESENTATION): $(PRESENTATION:%.pdf=%.tex) $(OUTPUT_DIR)
	$(DOC) bash -c "(TEXINPUTS=.///: xelatex -output-directory=$(OUTPUT_DIR) $<)"


$(OUTPUT_DIR)/$(TARGET): $(TARGET:%.pdf=%.tex) $(OUTPUT_DIR) $(KNITR) #$(KNITR_BASE)
	$(DOC) bash -c "(TEXINPUTS=.///: xelatex -output-directory=$(OUTPUT_DIR) $<) && (TEXINPUTS=.///: biber --output_directory=$(OUTPUT_DIR) $(<:%.tex=%)) && (TEXINPUTS=.///: xelatex -output-directory=$(OUTPUT_DIR) $<) && (TEXINPUTS=.///: xelatex -output-directory=$(OUTPUT_DIR) $<)"
$(OUTPUT_DIR):
	$(DOC) mkdir -p $(OUTPUT_DIR)

%.tex: %.Rnw $(KNITR)
	$(DOC) Rscript --vanilla -e "$(KNIT_COMMAND) knit('$(<F)',output='$(@F)')"

figures/%.tex: figures/%.Rnw 
	$(DOC) Rscript --vanilla -e "$(KNIT_COMMAND) knit('$(@D)/$(<F)', output = '$(@D)/$(@F)')"


superclean: clean
	$(DOC) find knitr -iname "*.tex" -delete
	$(DOC) find knitr -iname "*.pdf" -delete
	$(DOC) find knitr -iname "*.png" -delete
	$(DOC) find knitr -iname "images" -delete

clean:
	$(DOC) rm -rf  $(OUTPUT_DIR) $(TARGET) 


