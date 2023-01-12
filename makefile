all: process 

process: main_process.cpp tegbp.cpp io.cpp
	g++ -o process_gcc main_process.cpp tegbp.cpp plane_fitting.cpp io.cpp -fopenmp -Ofast
	# icpx -o process_intel  main_process.cpp tegbp.cpp plane_fitting.cpp io.cpp -qopenmp -Ofast

test: process 
	./process 32

clean:
	rm process_gcc 
	rm process_intel 
	# rm ./result/*
