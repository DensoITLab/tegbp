all: process 

process: main_process.cpp tegbp.cpp io.cpp
	g++ -o process main_process.cpp tegbp.cpp io.cpp -fopenmp -Ofast
	icpx -o process_i  main_process.cpp tegbp.cpp io.cpp -qopenmp -Ofast 
	# icc -o process -lm main_process.cpp tegbp.cpp io.cpp -fopenmp -Ofast
	#  icc -o process -lm main_process.cpp tegbp.cpp io.cpp -diag-disable=10441
	# g++ -o process -lm main_process.cpp tegbp.cpp io.cpp -fopenmp -Ofast -fexcess-precision=fast -ffast-math -mssse3 -mfpmath=sse

test: process 
	./process 32

clean:
	rm process 
	rm process_i 
	# rm ./result/*
