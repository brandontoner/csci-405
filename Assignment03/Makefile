TEX_FILES := $(wildcard *.tex)
PDF_FILES := $(TEX_FILES:.tex=.pdf)

ALL: $(PDF_FILES)

%.pdf: %.tex
	pdflatex $< < /dev/null

clean:
	rm -f *.pdf *.aux *.log

