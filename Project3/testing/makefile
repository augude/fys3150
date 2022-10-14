all: testPenningSetup testOneParticleFE testOneParticleRK4 testDoubleSetupWithInternal testDoubleSetupWithoutInternal compareStepsize

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
	for number in 4000 8000 16000 32000 ; do \
		g++ main.cpp src/tests.cpp ../../Project2/src/utils.cpp ../src/PenningTrap.cpp ../src/Particle.cpp -I include/ -larmadillo -o tests.exe  && ./tests.exe "CompareStepsize" $$number > $$number.txt ; \
	done