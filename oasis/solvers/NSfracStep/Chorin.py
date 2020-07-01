"""This is a simplest possible naive implementation of the Chorin solver.

The idea is that this solver can be quickly modified and tested for
alternative implementations. In the end it can be used to validate
the implementations of the more complex optimized solvers.

"""
from dolfin import *
import numpy as np
from scipy.sparse.linalg import splu
from scipy import sparse
from ..NSfracStep import *
from ..NSfracStep import __all__

__all__ += ["max_iter", "iters_on_first_timestep"]

# Chorin is noniterative
max_iter = 1
iters_on_first_timestep = 1


def setup(u, q_, q_1, uc_comp, u_components, dt, v, U_AB, u_1,
          nu, p_, dp_, mesh, f, fs, q, p, u_, Schmidt,
          scalar_components, **NS_namespace):
    """Set up all equations to be solved."""
    # Implicit Crank Nicholson velocity at t - dt/2
    U_CN = dict((ui, 0.5 * (u + q_1[ui])) for ui in uc_comp)

    F = {}
    Fu = {}
    for i, ui in enumerate(u_components):
        # Tentative velocity step
        F[ui] = ((1. / dt) * inner(u - q_1[ui], v) * dx
                 + inner(dot(U_AB, nabla_grad(U_CN[ui])), v) * dx
                 + nu * inner(grad(U_CN[ui]), grad(v)) * dx - inner(f[i], v) * dx)

        # Velocity update
        Fu[ui] = (inner(u, v) * dx - inner(q_[ui], v) *
                  dx + dt * inner(p_.dx(i), v) * dx)

    # Pressure solve
    Fp = inner(grad(q), grad(p)) * dx + (1. / dt) * div(u_) * q * dx

    # Scalar with SUPG
    h = CellDiameter(mesh)
    vw = v + h * inner(grad(v), U_AB)
    n = FacetNormal(mesh)
    for ci in scalar_components:
        F[ci] = ((1. / dt) * inner(u - q_1[ci], vw) * dx
                 + inner(dot(grad(U_CN[ci]), U_AB), vw) * dx
                 + nu / Schmidt[ci] * inner(grad(U_CN[ci]),
                            grad(vw)) * dx - inner(fs[ci], vw) * dx)

    return dict(F=F, Fu=Fu, Fp=Fp)


def velocity_tentative_solve(ui, F, q_, bcs, x_, b_tmp, udiff,
                             suppD, hm, Da, Db, dspac, rho, u_, dt, lu,
                             cor_x_Eul, cor_y_Eul, cor_z_Eul, **NS_namespace):
    """Linear algebra solve of tentative velocity component."""
    b_tmp[ui][:] = x_[ui]
    A, L = system(F[ui])
    solve(A == L, q_[ui], bcs[ui])
    udiff[0] += norm(b_tmp[ui] - x_[ui])
    
    ux = np.expand_dims(u_[0].vector().get_local(), axis=1)
    uy = np.expand_dims(u_[1].vector().get_local(), axis=1)
    uz = np.expand_dims(u_[2].vector().get_local(), axis=1)
    
    # Compute dU     
    Bibm_x = -(hm**3)*Db.dot(ux[suppD])
    Bibm_y = -(hm**3)*Db.dot(uy[suppD])
    Bibm_z = -(hm**3)*Db.dot(uz[suppD])
    
    dUx = lu.solve(Bibm_x)
    dUy = lu.solve(Bibm_y)
    dUz = lu.solve(Bibm_z)

    correction_x = dspac*Da.dot(dUx)
    correction_y = dspac*Da.dot(dUy)
    correction_z = dspac*Da.dot(dUz)
    
    cor_x_Eul[suppD] = correction_x
    cor_y_Eul[suppD] = correction_y
    cor_z_Eul[suppD] = correction_z
    
    # Updated velocity after IBM
    uc = ux + (dt/rho)*cor_x_Eul
    vc = uy + (dt/rho)*cor_y_Eul
    wc = uz + (dt/rho)*cor_z_Eul
    
    u_[0].vector()[:] = uc.ravel()
    u_[1].vector()[:] = vc.ravel()
    u_[2].vector()[:] = wc.ravel()
    
    for ui in u_components:
        [bc.apply(x_[ui]) for bc in bcs[ui]]


def pressure_solve(Fp, p_, bcs, **NS_namespace):
    """Solve pressure equation."""
    solve(lhs(Fp) == rhs(Fp), p_, bcs['p'])
    if bcs['p'] == []:
        normalize(p_.vector())


def velocity_update(u_components, q_, bcs, Fu, **NS_namespace):
    """Update the velocity after finishing pressure velocity iterations."""
    for ui in u_components:
        solve(lhs(Fu[ui]) == rhs(Fu[ui]), q_[ui], bcs[ui])


def scalar_solve(ci, F, q_, bcs, **NS_namespace):
    """Solve scalar equation."""
    solve(lhs(F[ci]) == rhs(F[ci]), q_[ci], bcs[ci])
