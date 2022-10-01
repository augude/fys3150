all: build run

build:
	g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe 

run:
	./tests.exe > ElectricField.txt