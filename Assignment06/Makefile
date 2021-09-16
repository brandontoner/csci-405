CXX       := g++
FLAGS     := -Wall -Wextra -pedantic -Werror -g -std=c++11
TEX_FILES := $(wildcard *.tex)
PDF_FILES := $(TEX_FILES:.tex=.pdf)
CPPLINT   := cpplint.py --filter=-whitespace/indent,-whitespace/semicolon,-runtime/threadsafe_fn,-runtime/printf,-legal/copyright

ALL: $(PDF_FILES) a6.out

a6.out: main.cpp
	dos2unix $^
	$(CPPLINT) $^
	$(CXX) $(FLAGS) -o $@ $^

clean:
	rm -f *.out *.pdf *.aux *.log

%.pdf: %.tex
	pdflatex $<

