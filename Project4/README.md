## Project 4 - The Ising model 
<p align = "center">
<img src = "https://user-images.githubusercontent.com/43797199/202687893-6a2c4188-31c4-4834-9eb0-c314ba46a022.png">
</p>
<p align = "center">
</b>In a remarkable stroke of luck we came across this microstate in one of our simulations, and thought we would include it as Ising on the cake.</b>
</p>


### Run code
To build and run the code for Project 4, write 

``make all``
in the terminal.
This will call commands stored in ``makefile``.

### Code structure
The header file ``lattice.hpp`` and the source file ``lattice.cpp`` contain a class for the Ising model. 
The header file ``mcmc.hpp`` and the source file ``mcmc.cpp`` contain the implementation of the Markov chain Monte Carlo algorithm.
Header files can be found in ``include`` and source files can be found in ``src``.
The main-function is build up with a series of functions that are each executed with a command-line argument. 
The folder ``plotting`` has notebooks for each problem that only serve plotting purposes. 
