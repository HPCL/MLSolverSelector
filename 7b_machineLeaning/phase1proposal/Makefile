# LaTeX M akefile for dvi, ps, and pdf file creation.
# By Jeffrey Humpherys
# Written April 05, 2004
# Revised January 13, 2005
# Thanks Bjorn and Boris

# Edited heavily by Gerald Sabin April 02, 2005

#
# Usage:
# make          # make dvi, ps, and pdf
# make dvi      # make dvi
# make ps       # make ps (and dvi)
# make pdf      # make pdf
#

MAIN            = main
SOURCES         = $(wildcard ./*.tex ./*/*.tex ./*.bib)
EPSFIGURES      = $(patsubst %.fig,%.ps,$(wildcard ./graphs/*.fig))

all : support Complete.pdf
summary: Summary.pdf
narrative: Narrative.pdf
company: Company.pdf
commercialization: Commercialization.pdf
budget: Budget.pdf
bio: Sabin_cv.pdf Robinson_cv.pdf Oneill_cv.pdf
complete: Complete.pdf
levelofeffort: LevelOfEffort.pdf

#Commercialization.pdf:
#	cp ./commercialization/commercial.pdf ./Commercialization.pdf

SBC.pdf:
	cp ./company/SBC.pdf ./SBC.pdf

Budget.pdf:	

Complete.pdf : Summary.pdf Narrative.pdf Sabin_cv.pdf Sabin_support.pdf ONeill_cv.pdf ONeill_support.pdf Robinson_cv.pdf Robinson_support.pdf Yallamelli_cv.pdf Yallamelli_support.pdf Budget.pdf LevelOfEffort.pdf Commercialization.pdf RNETHistory.pdf SBC.pdf Dataplan.pdf
	pdftk Summary.pdf Narrative.pdf Sabin_cv.pdf Sabin_support.pdf ONeill_cv.pdf ONeill_support.pdf Robinson_cv.pdf Robinson_support.pdf Budget.pdf LevelOfEffort.pdf Commercialization.pdf RNETHistory.pdf SBC.pdf Dataplan.pdf cat output Complete.pdf


ONeill_cv.pdf:
	cp ./bio/ben_cv_dec2017.pdf ONeill_cv.pdf

Robinson_cv.pdf:
	cp ./bio/jamar_cv.pdf Robinson_cv.pdf

Yallamelli_cv.pdf:
	cp ./bio/Pavan_CV.pdf Yallamelli_cv.pdf

Billings_cv.pdf:
	cp ./bio/billings_cv.pdf Billings_cv.pdf 

LevelOfEffort.pdf : 
	cp ./budget/level_of_effort.pdf ./LevelOfEffort.pdf

RNETHistory.pdf : Company.pdf
	pdftk Company.pdf ./company/RNET_Commercialization_History.pdf cat output RNETHistory.pdf
	rm Company.pdf

support : 
	pdftk ./bio/overall_support.pdf ./bio/robinson_support.pdf cat output ./Robinson_support.pdf
	pdftk ./bio/overall_support.pdf ./bio/oneill_support.pdf cat output ./Oneill_support.pdf
	pdftk ./bio/overall_support.pdf ./bio/sabin_support.pdf cat output ./Sabin_support.pdf
	pdftk ./bio/overall_support.pdf ./bio/yallamelli_support.pdf cat output ./Yallamelli_support.pdf

%.dvi : ${SOURCES} ${EPSFIGURES}
	./my_latex ${*}
	-bibtex ${*}
	./my_latex ${*}
	./my_latex ${*}
	@while ( grep "Rerun to get cross-references"                   \
		${*}.log .log> /dev/null ); do                   \
		echo '** Re-running LaTeX **';                          \
		./my_latex ${*};                                          \
	done

%.pdf : ${SOURCES} ${EPSFIGURES}
	./my_pdflatex ${*}
	-bibtex ${*}
	./my_pdflatex ${*}
	./my_pdflatex ${*}

#%.ps : %.dvi
#	dvips ${<} -o ${@}

%.rtf : %.dvi
	latex2rtf ${*}
	#echo "skipping rtf"

sched.bib:

clean:
	rm -f ./*.log
	rm -f ./*.aux
	rm -f ./*.out
	rm -f ./*.bbl
	rm -f ./*.blg
	rm -f ./*~
	rm -f ./*.bak
	rm -f ./\#*\#
	rm -f ./*/*~
	rm -f ./*/*.bak
	rm -f ./*/\#*\#
	rm -f ./*/*/*~
	rm -f ./*/*/*.bak
	rm -f ./*/*/\#*\#

veryclean: clean
	rm -f *.ps
	rm -f *.dvi
	rm -f *.rtf

pdfclean: veryclean
	rm *.pdf
