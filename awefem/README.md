# awefem
Time domain wave equation solver

Description
-----------
The script [`awefem.py`](awefem.py) solves the constant velocity scalar wave equation in an arbitrary number of dimensions. It injects a point source with a time-dependent source time function. Templates for a Ricker and sine waves are used. Dirichlet boundary conditions set to zero are used.

Disk
--------
The unit disk is discretized using a resolution of 50 elements. Source time function is a 40Hz sine function at the center of the disk.

Wireframe | 2D | 3D
----------|----|----
![Time evolution of the wave equation](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/circle/circle.gif)|![Time evolution of the wave equation](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/circle/circle_flat.gif)|![Time evolution of the wave equation](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/circle/circle_wire.png)



