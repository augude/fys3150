all: problem4 problem 5

server4:
	g++ -std=c++11 main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && ./main.exe "validation2x2"

server5:
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "burnIn"

server6:
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "histograms"

server7:
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe 
	for num in 1 2 3 4 ; do \
		export OMP_NUM_THREADS=$$num && ./main.exe "par" >> "problem7.txt"; \
	done

server8a: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && export OMP_NUM_THREADS=32 && ./main.exe "problem8" 40

server8b: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && export OMP_NUM_THREADS=32 && ./main.exe "problem8" 60

server8c: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && export OMP_NUM_THREADS=32 && ./main.exe "problem8" 80

server8d:
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && export OMP_NUM_THREADS=32 && ./main.exe "problem8" 100

server8bonus:
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -fopenmp -o main.exe && export OMP_NUM_THREADS=32 && ./main.exe "problem8" 1000


problem4: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && ./main.exe "validation2x2"

problem5: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "burnIn"

problem6: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "histograms"

problem7:
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe 
	for num in 1 2 3 4 ; do \
		export OMP_NUM_THREADS=$$num && ./main.exe "par" >> "problem7.txt"; \
	done

problem8a: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=8 && ./main.exe "problem8" 40

problem8b: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=8 && ./main.exe "problem8" 60

problem8c: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=8 && ./main.exe "problem8" 80

problem8d: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=8 && ./main.exe "problem8" 100

problem8: problem8a problem8b problem8c problem8d

problem8aZ: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "problem8zoom" 40

problem8bZ: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "problem8zoom" 60

problem8cZ: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "problem8zoom" 80

problem8dZ: 
	g++ main.cpp src/lattice.cpp src/mcmc.cpp src/utils.cpp -O2 -I \include -larmadillo -fopenmp -o main.exe && export OMP_NUM_THREADS=4 && ./main.exe "problem8zoom" 100

problem8Z: problem8aZ problem8bZ problem8cZ problem8dZ
