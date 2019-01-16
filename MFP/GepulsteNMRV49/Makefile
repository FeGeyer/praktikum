all: build/dokument.pdf

TeXOptions = -lualatex \
				 -interaction=nonstopmode \
				 -halt-on-error \
				 -output-directory=build

build/dokument.pdf: FORCE | build
	mkdir -p Auswertung/Tabellen/
	mkdir -p Auswertung/Plots/
	python 	auswertung.py
	latexmk $(TexOptions) dokument.tex

light: dokument.tex
	lualatex -interaction=nonstopmode \
	-halt-on-error \
	-output-directory=build \
	dokument.tex

FORCE:

build:
	mkdir -p build/

clean:
	rm -rf build
