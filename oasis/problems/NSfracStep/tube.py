from dolfin import *
from mshr import *
from ..NSfracStep import *

from scipy.sparse.linalg import splu
from scipy import sparse
import numpy as np

try:
    from matplotlib import pyplot as plt
except:
    pass

# Override some problem specific parameters
def problem_parameters(NS_parameters, commandline_kwargs, NS_expressions, **NS_namespace):
    NS_parameters.update(
        nu=0.00345/1056,
        rho=1056,
        T=1.0,
        dt=0.0005,
        max_iter=1,                 # Number of inner pressure velocity iterations on timestep
        iters_on_first_timestep=3,  # Number of iterations on first timestep
        folder="Immersed boundary tube",
        plot_interval=500,
        save_step=20,
        checkpoint=5000,
        print_intermediate_info=10,
        compute_error=1,
        use_krylov_solvers=True,
        velocity_degree=1,
        pressure_degree=1,
        krylov_report=False)
    
    if "velocity_degree" in commandline_kwargs.keys():
        degree = commandline_kwargs["velocity_degree"]
    else:
        degree = NS_parameters["velocity_degree"]

    NS_expressions.update(dict(u_in = Expression('0.122625*(1 - (x[0]*x[0] + x[2]*x[2])/(0.004*0.004))', degree=degree)))

    NS_parameters['krylov_solvers'] = {'monitor_convergence': False,
                                       'report': False,
                                       'relative_tolerance': 1e-12,
                                       'absolute_tolerance': 1e-12}

# Create a mesh
def mesh(**params):
    #Eulerian_mesh 
    m = Mesh('/home/usr/mesh_ref_250.xml')
    return m

def create_bcs(mesh, V, Q, sys_comp, u_in, **NS_namespace):
    
    # Define boundaries
    center = Point(0.0, 0.0)
    radius = 0.004 # Cylinder radious
    
    #h=0.00025 m
    Lxmax = 0.032
    Lymin = 0.0

    # No-slip boundary (walls)
    class Walls(SubDomain):
        def inside(self, x, boundary):
            return boundary 

    # Inflow boundary
    class Inflow(SubDomain):
        def inside(self, x, boundary):
            oncirclehaut = np.hypot(x[0] - center[0], x[2] - center[1]) 
            return oncirclehaut <= radius and near(x[1], Lymin) and boundary
            
    # Outflow boundary
    class Outflow(SubDomain):
        def inside(self, x, boundary):
            return near(x[0], Lxmax) and boundary
    
    # Create mesh functions over the cell facets
    sub_domains = MeshFunction("size_t", mesh, mesh.topology().dim() - 1)

    # Mark all facets as sub domain 3
    sub_domains.set_all(0)

    # Mark no-slip facets as sub domain 0
    noslip = Walls()
    noslip.mark(sub_domains, 1)

    # Mark inflow as sub domain 1
    inflow = Inflow()
    inflow.mark(sub_domains, 2)

    # Mark outflow as sub domain 2, 
    outflow = Outflow()
    outflow.mark(sub_domains, 3)    

    bcs = dict((ui, []) for ui in sys_comp)
    bc0 = DirichletBC(V, 0, sub_domains, 1)
    bc1 = DirichletBC(V, u_in, sub_domains, 2)
    bc2 = DirichletBC(V, 0, sub_domains, 2)
    pc1 = DirichletBC(Q, 0.0, sub_domains, 3)
    bcs['u0'] = [bc0, bc2]
    bcs['u1'] = [bc0, bc1]
    bcs['u2'] = [bc0, bc2]
    bcs['p'] = [pc1]
    
    return bcs
                
def Eul_nodes(V, mesh, **NS_namespace):
    
    Eulerian = V.tabulate_dof_coordinates()    
    hm = np.abs(Eulerian[1,0] - Eulerian[0,0])
    hm = hm / 16
    return dict(Eulerian=Eulerian, hm=hm)

def ibmat( **NS_namespace):
    
    """Import the (sparse) immersed boundary matrices."""
        
    print("Importing the IB matrices")
    
    dspac = np.load('/home/usr/dspac.npy')
    dspac = np.mean(dspac)

    suppD = np.load('/home/usr/suppD.npy')
    
    Da = sparse.load_npz('/home/usr/Meshes/Da.npz')
    Db = sparse.load_npz('/home/usr/Meshes/Db.npz')

    return dict(dspac=dspac, Da=Da, Db=Db, suppD=suppD)

def Aib(suppD, Da, Db, rho, dt, hm, dspac, Eulerian, **NS_namespace):
    
    Aibm = dt*(1/rho)*(hm**3)*((Db).dot(Da*dspac))
    lu = splu(Aibm)

    cor_x_Eul = np.zeros((Eulerian.shape[0],1))
    cor_y_Eul = np.zeros((Eulerian.shape[0],1))
    cor_z_Eul = np.zeros((Eulerian.shape[0],1))

    da = dict(lu=lu, Aibm=Aibm, cor_x_Eul=cor_x_Eul, cor_y_Eul=cor_y_Eul, cor_z_Eul=cor_z_Eul)    
    
    return da

def initialize(x_1, x_2, bcs, **NS_namespace):
    for ui in x_1:
        [bc.apply(x_1[ui]) for bc in bcs[ui]]
    for ui in x_2:
        [bc.apply(x_2[ui]) for bc in bcs[ui]]

def pre_solve_hook(mesh, velocity_degree, u_,
                   AssignedVectorFunction, **NS_namespace):
    return dict(uv=AssignedVectorFunction(u_))

def temporal_hook(tstep, uv, u_, plot_interval, **NS_namespace):
    
    if tstep % plot_interval == 0:
        
        # Assign solution to vector
        uv()

        
