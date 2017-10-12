# IsingModel
This project focuses on general Ising model. 
- Designed for general spins.
- Set spin, calculate energy and then flip. Example is given for square lattice spin half, nearest neighbourhood and single flip
- Above three functions should be overrode for other models



For openmp
- do remember to compile with -fopenmp
- TDM-GCC by default does not include openmp! Remember to check
- If the Compiler Flag Options is empty (Code::Blocks), delete the xml in AppData\C::B\share\....\....\compiler
- Add a Flag with -fopenmp, and add libgomp*.dll as link library