from dolfin import UnitSquareMesh
from numpy import cos, pi


# Create a mesh
def mesh(Nx=160, Ny=160, **params):
    m = UnitSquareMesh(Nx, Ny)
    x = m.coordinates()
    #x[:] = (x - 0.5) * 2
    #x[:] = 0.5*(cos(pi*(x-1.) / 2.) + 1.)
    return m


noslip = "std::abs(x[0]*x[1]*(1-x[0]))<1e-8"
top = "std::abs(x[1]-1) < 1e-8"
bottom = "std::abs(x[1]) < 1e-8"
