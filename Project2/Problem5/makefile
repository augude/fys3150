all: build run

build:
	g++ main.cpp src/countSimilarity.cpp ../Problem4/src/jacobi.cpp ../src/utils.cpp ../Problem2/src/setupTridiag.cpp ../Problem3/src/maxSymMatrix.cpp -I include/ -larmadillo -o countSimilarities.exe
run:
	./countSimilarities.exe