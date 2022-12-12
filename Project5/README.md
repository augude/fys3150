## <p align = "center">Project 5 - Quantum double-slit </p>

### Run code
To build and run the code for Project 5, write ``make all`` in the terminal. This will call commands stored in ``makefile``.
The python script ``create_input.py`` in the folder ``input`` will create the matrix holding the values for the potential and save the output to a ``dat``-file.
The output of the simulation is stored in ``bin``-files and can be visualized using the functions in the ``plotting_utils.py`` in the folder ``plotting``.


### Code structure
The header file ``Solver.hpp`` and the source file ``Solver.cpp`` contain a class for the Crank-Nicholson solver. 
Header files can be found in ``include`` and source files can be found in ``src``.
The folder ``plotting`` has a file with stored function for visualization and a notebook that serve plotting purposes.
Some animations of the time-evolution can be found in the folder ``output``.
