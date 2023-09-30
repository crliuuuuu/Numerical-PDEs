run:	TI.h test.cpp test.json
	g++ -o test test.cpp -ljsoncpp -std=c++17 -O3 -Ofast
	./test
story:	./report_src/report.tex
	xelatex ./report_src/report.tex
	rm -f report.aux
	rm -f report.bbl
	rm -f report.log
	rm -f report.blg
debug:	TI.h test.cpp test.json
	g++ -o test test.cpp -ljsoncpp -std=c++17 -g
clean:
	rm -f test report.pdf
