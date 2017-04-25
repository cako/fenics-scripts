#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

from __future__ import print_function, division
from dolfin import *
from mshr import Circle, generate_mesh
import numpy as np
import os, sys, time

# Set log level
set_log_active(False)

def ricker_source(t, f=40):
    t -= 5*np.sqrt(3/2) / (f*np.pi)
    return (1 - 2*(np.pi*f*t)**2) * np.exp( -(np.pi*f*t)**2 )

def sine_source(t, f=40):
    return np.sin(2*np.pi*f*t)

def awefem(mesh, t, source_time_function = None, source_loc = None,
           filename=None, verb=False):
    start_time = time.time()

    # Function space
    V = FunctionSpace(mesh, "Lagrange", 1)

    # Boundary condition
    # TODO: Absorbing boundary conditions
    bc = DirichletBC(V, Constant(0), 'on_boundary')

    # Trial and test functions
    u = TrialFunction(V)
    v = TestFunction(V)

    # Discretization
    # ∂tt u = c²Δu
    # (uN+1 - 2uN + uN-1)/dt² = c²Δu
    # ∫ (uN+1 - 2uN + uN-1)v dx - (dt*c)²Δu v dx
    # ∫ (uN+1 - 2uN + uN-1)v dx + (dt*c)²∇u·∇v dx
    # ∫ (uN+1 - 2uN + uN-1)v dx + (dt*c)²∇(uN+1 + 2uN + uN-1)·∇v dx
    c = 6
    dt, nt = t[1]-t[0], len(t)
    u0 = Function(V, name = "u0") # u0 = uN-1
    u1 = Function(V, name = "u1") # u1 = uN1

    # Variational formulation
    F = (u-2*u1+u0)*v*dx + (dt*c)**2*dot(grad(u+2*u1+u0)/4, grad(v))*dx
    a, L = lhs(F), rhs(F)

    # Solver
    A, b  = assemble_system(a, L)
    solver = LUSolver(A)
    solver.parameters['reuse_factorization'] = True
    bc.apply(A,b)

    # Solution
    u = Function(V, name = "u") # uN+1

    # Source
    if source_time_function is None:
        source_time_function = ricker_source
    if source_loc is None:
        mesh_center = np.mean(mesh.coordinates(), axis=0)
        source_loc = Point(mesh_center)

    # Aux
    if filename is not None:
        vtkfile = File(filename)
    bk = 2*len(str(nt))+1

    # Time stepping
    for i, t_ in enumerate(t[1:]):
        sys.stdout.write('\010'*bk)
        sys.stdout.flush()
        print('%d/%d' % (i+2, nt), end='')

        b  = assemble(L)
        delta = PointSource(V, source_loc, source_time_function(t_)*dt**2)
        delta.apply(b)
        solver.solve(u.vector(), b)

        if filename is not None:
            vtkfile << (u, t_)

        u0.assign(u1)
        u1.assign(u)
    print(time.strftime("\nElapsed time %H:%M:%S",
                        time.gmtime(time.time()-start_time)))

    return u

if __name__ == "__main__":
    ot, dt, nt = 0., 5e-4, 1001
    t = ot + np.arange(nt)*dt

    print('Computing wavefields over dolfin mesh')
    mesh = Mesh("../meshes/dolfin_fine.xml.gz")
    u = awefem(mesh, t, source_loc = Point(0.8, 0.8),
            filename='dolfin/data/dolfin.pvd')

    print('Computing wavefields over unit square')
    mesh = UnitSquareMesh(100, 100)
    u = awefem(mesh, t, filename='square/data/square.pvd')

    print('Computing wavefields over unit circle')
    domain = Circle(Point(0., 0.), 1)
    mesh = generate_mesh(domain, 50)
    u = awefem(mesh, t, source_time_function = sine_source,
               filename='circle/data/circle.pvd')
