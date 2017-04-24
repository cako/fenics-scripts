# awefem
Time domain wave equation solver

Description
-----------
The script [`awefem.py`](awefem.py) solves the constant velocity scalar wave equation in an arbitrary number of dimensions. It injects a point source with a time-dependent source time function. Templates for a Ricker and sine waves are used. Dirichlet boundary conditions set to zero are used.

Square
------
Unit square mesh with a Ricker wavelet in the center.

![Time evolution of the wave equation](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/circle/circle.gif)
