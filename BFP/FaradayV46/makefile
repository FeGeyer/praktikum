ifeq (,$(shell sh -c 'cygpath --version 2> /dev/null'))
  # Unix
  pwd := $$(pwd)
  translate = $1
else
  # Windows mit MSys2/Cygwin
  pwd := $$(cygpath -m "$$(pwd)")
  translate = $(shell echo '$1' | sed 's/:/;/g')
endif

all: build/dokument.pdf

# hier Python-Skripte:
# build/plot.pdf: auswertung.py matplotlibrc header-matplotlib.tex | build
# 	TEXINPUTS="$(call translate,$(pwd):)" python auswertung.py

# hier weitere Abhängigkeiten für build/main.pdf deklarieren:
# build/dokument.pdf: build/plot.pdf

build/dokument.pdf: FORCE | build
	  TEXINPUTS="$(call translate,build:)" \
	  BIBINPUTS=build: \
	  max_print_line=1048576 \
	latexmk \
	  --lualatex \
	  --output-directory=build \
	  --interaction=nonstopmode \
	  --halt-on-error \
	dokument.tex

build:
	mkdir -p build

clean:
	rm -rf build

FORCE:

.PHONY: all clean
