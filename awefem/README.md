# awefem
Time domain wave equation solver

Description
-----------
The script [`awefem.py`](awefem.py) solves the constant velocity scalar wave equation in an arbitrary number of dimensions. It injects a point source with a time-dependent source time function. Templates for a Ricker and sine waves are used. Dirichlet boundary conditions set to zero are used.

### Disk
The unit disk is discretized using a resolution of 50 elements. Source time function is a 40Hz sine function at the center of the disk.

Wireframe | 2D | 3D
----------|----|----
![Wireframe of domain](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/circle/circle_wire.png)|![Time evolution of the wave equation 2D](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/circle/circle_flat.gif)|![Time evolution of the wave equation in 3D](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/circle/circle.gif)

### Dolphin
The unit square with a dolpin shape remove. The mesh is unstructured as can be seen below. Source time function is a 40Hz Ricker wavelet near the top right edge.

Wireframe | 2D | 3D
----------|----|----
![Wireframe of domain](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/dolfin/dolfin_wire.png)|![Time evolution of the wave equation 2D](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/dolfin/dolfin_flat.gif)|![Time evolution of the wave equation in 3D](https://raw.githubusercontent.com/cako/fenics-scripts/master/awefem/dolfin/dolfin.gif)


