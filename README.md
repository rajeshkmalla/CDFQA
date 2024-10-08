# Feedback-based Quantum Algorithm Inspired by Counterdiabatic Driving
The codes for classical simulations are written in Matlab and then converted to Python. We include both Matlab codes as well as Python codes for convenience.

## Matlab codes

The *HamNN* function constructs the Hamiltonian for a 1D quantum spin chain with N sites and nearest-neighbor interactions.

The *HamOnsite* function constructs the Hamiltonian for a 1D quantum spin chain with on-site interactions (i.e., individual spin terms at each site). The spin operator (spin1) is applied at every site, and the matrix size scales with N (the number of sites).

The *fastExpm* function efficiently computes the matrix exponential, designed to work with both sparse and full matrices, and even on GPU if available. It uses a combination of scaling, Taylor series expansion, and squaring for performance, while maintaining sparsity when possible.
This code has been originally developed by Ilya Kuprov (http://spindynamics.org/) and has been adapted by F. Mentink-Vigier (fmentink@magnet.fsu.edu)
% If you use this code, please cite 
%  - H. J. Hogben, M. Krzystyniak, G. T. P. Charnock, P. J. Hore and I. Kuprov, Spinach – A software library for simulation of spin dynamics in large spin systems, J. Magn. Reson., 2011, 208, 179–194.
%  - I. Kuprov, Diagonalization-free implementation of spin relaxation theory for large spin systems., J. Magn. Reson., 2011, 209, 31–38.

** All the Matlab functions and the file to generate energy vs Cir depth have been converted to Python. **

## Python codes for executing on quantum computers

The Run_on_quantum_computer.ipynb file includes constructing the feedback-based circuits and the procedure for running the circuits iteratively on quantum computers. The codes are mostly using the Qiskit package (https://github.com/qiskit) with version 1.2.0. The file itself only contains the methods and results for running on a local simulator. 

To execute the circuits with real quantum computers, change the variable "backend" accordingly.

To add error mitigation methods, modify the flags under the variable "options" accordingly.


## Python codes for Pennylane
  
Includes example codes for simulation on Pennylane
