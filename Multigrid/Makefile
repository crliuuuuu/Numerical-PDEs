run:	MC.h test.cpp test.json
	g++ -o test test.cpp -ljsoncpp -std=c++11 -O3
	./test
debug:	MC.h test.cpp test.json
	g++ -o test test.cpp -ljsoncpp -std=c++11 -g
clean:
	rm -f test 
