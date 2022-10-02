all: testPenningSetup testOneParticleFE testOneParticleRK4 testDoubleSetupWithInternal testDoubleSetupWithoutInternal

testPenningSetup: 
	g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe  && ./tests.exe "PenningSetup" > ElectricField.txt	

testOneParticleFE: 
	g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe  && ./tests.exe "OneParticleFE" > testOneParticleFE.txt	

testOneParticleRK4: 
	g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe  && ./tests.exe "OneParticleRK4" > testOneParticleRK4.txt	

testDoubleSetupWithInternal: 
	g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe  && ./tests.exe "DoubleSetupWithInternal" > testDoubleSetupWithInternal.txt	

testDoubleSetupWithoutInternal: 
	g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe  && ./tests.exe "DoubleSetupWithoutInternal" > testDoubleSetupWithoutInternal.txt	

compareStepsize:
	for number in 1e-1 1e-2 1e-3 1e-4 1e-5 ; do \
		g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe  && ./tests.exe "CompareStepsize" $$number > $$number.txt ; \
	done