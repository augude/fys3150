all: problem4 problem 5

problem4: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && ./main.exe "validation2x2"

problem5: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "burnIn"

problem6: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=5 && ./main.exe "histograms"

problem7:
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe 
	for num in 1 2 3 4 ; do \
		export OMP_NUM_THREADS=$$num && ./main.exe "par" >> "problem7.txt"; \
	done

problem8a: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=10 && ./main.exe "problem8" 40

problem8b: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=10 && ./main.exe "problem8" 60

problem8c: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=10 && ./main.exe "problem8" 80

problem8d: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=10 && ./main.exe "problem8" 100
