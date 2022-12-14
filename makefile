all: process 

process: main_process.cpp tegbp.cpp io.cpp
	g++ -o process -lm main_process.cpp tegbp.cpp io.cpp -fopenmp -Ofast
	# g++ -o process -lm main_process.cpp tegbp.cpp io.cpp -fopenmp -Ofast -fexcess-precision=fast -ffast-math -mssse3 -mfpmath=sse

test: process 
	./process 32

clean:
	rm process 
