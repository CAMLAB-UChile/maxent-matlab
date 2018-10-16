# maxent-matlab
A MATLAB implementation of the maximum-entropy basis functions

# Author
<a href="https://github.com/aaortizb">Alejandro Ortiz-Bernardin</a>

# Instructions
The program is controlled by the main.m function. This is the only function that 
must be setup by the user. To execute the code, setup the problem parameters in
main.m (further instructions are given there) and run it.

When setting up main.m make sure that
  - size(x,2) = dim
  - length(x) = size(ncoord,2)
  - size(ncoord,1) = n
  - length(gamma) = n
  - length(ilambda) = dim
  
Anyway, an error is thrown when any of the previous equalities are not satisfied.

# License
This project is licensed under the GPL3 License. This program is free software; it can be redistributed or modified under the terms of the GNU General Public License 3 as published by the Free Software Foundation. 
