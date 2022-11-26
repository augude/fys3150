noslits:
	cd input \
	&& python3 create_input.py no 5e-3 0 \
	&& cd .. && g++ main.cpp src/Solver.cpp -O2 -I \include -larmadillo -o main.exe && ./main.exe 5e-3 2.5e-5 2e-3 0.25 0.05 200 0.5 0.05 0 no

doubleslits_long:
	cd input \
	&& python3 create_input.py double 5e-3 1e10 \
	&& cd .. && g++ main.cpp src/Solver.cpp -O2 -I \include -larmadillo -o main.exe && ./main.exe 5e-3 2.5e-5 8e-3 0.25 0.05 200 0.5 0.10 0 double

doubleslits_short:
	cd input \
	&& python3 create_input.py double 5e-3 1e10 \
	&& cd .. && g++ main.cpp src/Solver.cpp -O2 -I \include -larmadillo -o main.exe && ./main.exe 5e-3 2.5e-5 2e-3 0.25 0.05 200 0.5 0.20 0 double

singleslit:
	cd input \
	&& python3 create_input.py single 5e-3 1e10 \
	&& cd .. && g++ main.cpp src/Solver.cpp -O2 -I \include -larmadillo -o main.exe && ./main.exe 5e-3 2.5e-5 2e-3 0.25 0.05 200 0.5 0.20 0 single


tripleslits:
	cd input \
	&& python3 create_input.py triple 5e-3 1e10 \
	&& cd .. && g++ main.cpp src/Solver.cpp -O2 -I \include -larmadillo -o main.exe && ./main.exe 5e-3 2.5e-5 2e-3 0.25 0.05 200 0.5 0.20 0 triple


tunnel:
	cd input \
	&& python3 create_input.py tunnel 5e-3 4e4 \
	&& cd .. && g++ main.cpp src/Solver.cpp -O2 -I \include -larmadillo -o main.exe && ./main.exe 5e-3 2.5e-5 8e-3 0.25 0.05 200 0.5 0.20 0 tunnel