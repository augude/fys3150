all: problem4 problem 5

problem4: 
	g++ main.cpp src/utils.cpp src/lattice.cpp src/mcmc.cpp -O2 -I \include -larmadillo -o main.exe && ./main.exe "validation2x2"

problem5: 
	g++ main.cpp src/utils.cpp src/lattice.cpp src/mcmc.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "burnIn"