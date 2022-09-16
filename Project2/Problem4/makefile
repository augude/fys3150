all: build run

build:
	g++ main.cpp src/jacobi.cpp ../src/utils.cpp ../Problem3/src/maxSymMatrix.cpp ../Problem2/src/setupTridiag.cpp -I include/ -larmadillo -o jacobi.exe
run:
	./jacobi.exe