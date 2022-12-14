all: process 

process: main_process.cpp tegbp.cpp io.cpp
	g++ -o process -lm main_process.cpp tegbp.cpp io.cpp -fopenmp -O3

test: process 
	./process 32

clean:
	rm process 
