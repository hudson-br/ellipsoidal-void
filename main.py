from dolfin import *
from mshr import *    
import meshio
import gmsh
import matplotlib.pyplot as plt
import numpy as np
import os, sympy, shutil
import ufl
os.environ['KMP_DUPLICATE_LIB_OK']='True'

parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["quadrature_degree"] = 2
parameters["form_compiler"]["representation"]="uflacs"
parameters["use_petsc_signal_handler"] = True
info(parameters,True)


R = 1
R_void = 0.2

R_void2 = 3*R_void 
distance_between_voids = 0 * R_void


load0 = 1. # reference value for the loading (imposed displacement)
load_multipliers = load0*np.linspace(0,.1,101)

filename =  "domain"

savedir = "result/"
if not os.path.exists(savedir):
    os.makedirs(savedir)
if os.path.isdir(savedir):
    shutil.rmtree(savedir)   
    output_file = XDMFFile(savedir + "results.xdmf")
    output_file.parameters["functions_share_mesh"] = True
    output_file.parameters["flush_output"] = True





    
# Leaf directory 
directory = "last_iteration_damage"
    
# Parent Directories 
parent_dir = os.getcwd()
    
# Path 
path = os.path.join(parent_dir, directory) 
    
# Create the directory 

try:
    os.makedirs(path, exist_ok = True)
    print("Directory '%s' created successfully" % directory)
except OSError as error:
    print("Directory '%s' can not be created" % directory)




# Read mesh


xdmf_name = filename + ".xdmf"

# msh = meshio.read(xdmf_name.replace(".xdmf", ".msh"))
msh = meshio.read(filename+ ".msh")

meshio.write(
    xdmf_name,
    meshio.Mesh(points=msh.points[:,:2], cells={"triangle": msh.cells_dict["triangle"]}),
)
mmesh = meshio.read(xdmf_name)
mesh = Mesh()
with XDMFFile(xdmf_name) as mesh_file:
    mesh_file.read(mesh)





exterior_circle = CompiledSubDomain("near(pow(x[0],2) + pow(x[1],2), pow(%f,2), 1.e-2)"%R)

# hole1 = CompiledSubDomain("near(pow(x[0] + 0.2,2) + pow(x[1],2), pow(%f,2), 1.e-3)"%R )
hole1 = CompiledSubDomain("near(pow(x[0] - {R0},2) + pow(x[1],2), pow({R},2), 1.e-3)".format(R0 = -distance_between_voids/2, R = R_void) )
# hole2 = CompiledSubDomain("near(pow(x[0] - 0.2,2) + pow(x[1],2), pow(%f,2), 1.e-3)"%R )
hole2 = CompiledSubDomain("near(pow(x[0] - {R0},2) + pow(x[1],2), pow({R},2), 1.e-3)".format(R0 =  distance_between_voids/2, R = R_void) )
# hole3 = CompiledSubDomain("near(pow(x[0] - {R0},2) + pow(x[1],2), pow({R},2), 1.e-3)".format(R0 = R03, R = R3) )

boundaries = MeshFunction("size_t", mesh,1)
boundaries.set_all(0)
exterior_circle.mark(boundaries, 1)
hole1.mark(boundaries, 2)
hole2.mark(boundaries, 3)
# hole3.mark(boundaries, 7)



# bottom.mark(boundaries, 2)
ds = Measure("ds",subdomain_data=boundaries) 

ndim = mesh.topology().dim() # get number of space dimensions
zero_v = Constant((0.,)*ndim) # a ndim-dimensional zero vector
print(mesh.hmax())

cell_markers = MeshFunction("bool", mesh, 1)
cell_markers.set_all(False)


normals = FacetNormal(mesh)

####################################






E = Constant(2.0)
Gc = Constant(1.0)
ell = Constant(0.04)

mu = 0.5*E



def w(alpha):
    """Dissipated energy function as a function of the damage """
    return alpha

def a(alpha):
    """Stiffness modulation as a function of the damage """
    k_ell = Constant(1.e-4) # residual stiffness
    return (1-alpha)**2+k_ell


z = sympy.Symbol("z")
c_w = 4*sympy.integrate(sympy.sqrt(w(z)),(z,0,1))
print("c_w = ",c_w)
c_1w = sympy.integrate(sympy.sqrt(1/w(z)),(z,0,1))
print("c_1/w = ",c_1w)
tmp = 2*(sympy.diff(w(z),z)/sympy.diff(1/a(z),z)).subs({"z":0})
sigma_c = sympy.sqrt(tmp*Gc*E/(c_w*ell))
print("sigma_c = %2.3f"%sigma_c)
eps_c = float(sigma_c/E)
print("eps_c = %2.3f"%eps_c)


#####################################

# Create function space for 2D elasticity + Damage
P2 = FiniteElement("Lagrange", mesh.ufl_cell(), degree=1)

element = MixedElement([VectorElement(P2, dim=2),  P2])
V = FunctionSpace(mesh, element)
V_tensor = TensorFunctionSpace(mesh, "CG", 1, shape=(2, 2))
V_alpha = FunctionSpace(mesh, "P", 1)



# Define the function, test and trial fields
q, du, v = Function(V), TrialFunction(V), TestFunction(V)
u, p = split(q)
alpha, dalpha, beta = Function(V_alpha), TrialFunction(V_alpha), TestFunction(V_alpha)

dissipated_energy = Gc/float(c_w)*(w(alpha)/ell+ ell*dot(grad(alpha), grad(alpha)))*dx

# Weak form of elasticity problem
# E_u = derivative(total_energy,u,v)

###########################

# expression to move inner surface x component
# xmove = Expression("t*x[0]/(pow(pow(x[0],2) + pow(x[1],2), 0.5))", t = 0,  degree=2)
xmove = Expression("t*(x[0] - {R0}) /(pow(pow(x[0] - {R0},2) + pow(x[1],2), 0.5))".format(R0 = 0, R1 = R), t = 0,  degree=2)
# expression to move inner surface y component
ymove = Expression("t*x[1]/(pow(pow(x[0] - {R0},2) + pow(x[1],2), 0.5))".format(R0 = 0, R1= R), t = 0,  degree=2)


bc_ic_x = DirichletBC(V.sub(0).sub(0), xmove, boundaries, 1) # u0 move - Inside
bc_ic_y = DirichletBC(V.sub(0).sub(1), ymove, boundaries, 1) # u0 move - Inside


# Dirichlet boundary condition for a traction test boundary
bc_u = [bc_ic_x, bc_ic_y]
# bc_u =[]


bcalpha_1 = DirichletBC(V_alpha, 0.0, boundaries, 1)

bc_alpha = [bcalpha_1]



#####################################







# E_du = ufl.replace(E_u,{u:du})
# problem_u = LinearVariationalProblem(lhs(E_du), rhs(E_du), u, bc_u)
# solver_u = LinearVariationalSolver(problem_u)
# solver_u.parameters.update({"linear_solver" : "mumps"})




# Identity tensor
d = len(u) # Spatial dimension

I = ufl.variable(ufl.Identity(d))
x = SpatialCoordinate(mesh)

# Deformation gradient

F = ufl.variable(ufl.grad(u) + I)


# Right Cauchy-Green tensor
C = ufl.variable(F.T * F)
# Left Cauchy-Green tensor
B = ufl.variable(F * F.T)

# Invariants of deformation tensors
Ic = ufl.variable(ufl.tr(C))
J  = ufl.variable(ufl.det(F))


T = ufl.variable(-p*I + mu*(B-I)) # Cauchy stress
S = J*T*inv(F).T # 1st Piola-Kirchhoff stress

# Stored strain energy density (compressible neo-Hookean model)
# psi = ufl.variable((mu / 2.0) * (Ic - 2) - mu * ufl.ln(J) + (lmbda / 2.0) * (ufl.ln(J))**2)


psi = ufl.variable((mu / 2.0) * (Ic - 2) - p*(ufl.det(F) - 1)) # - mu * ufl.ln(J) + (lmbda / 2.0) * (ufl.ln(J))**2)


Pi = a(alpha) * psi * dx # + p*dot(u, n)*ds(2) + p*dot(u, n)*ds(3)



    
F_res = ufl.derivative(Pi, q, v)
J = ufl.derivative(F_res, q, du)

total_energy = Pi + dissipated_energy 
E_alpha = derivative(total_energy,alpha,beta)
E_alpha_alpha = derivative(E_alpha,alpha,dalpha)



problem_unl = NonlinearVariationalProblem(F_res, q, bc_u, J)
solver_u = NonlinearVariationalSolver(problem_unl)
prm = solver_u.parameters
prm["newton_solver"]["error_on_nonconvergence"] = True
prm["newton_solver"]["absolute_tolerance"] = 1e-3
prm["newton_solver"]["relative_tolerance"] = 1e-3
# solver.solve()

##########################################

class DamageProblem(OptimisationProblem):

    def f(self, x):
        """Function to be minimised"""
        alpha.vector()[:] = x
        return assemble(total_energy)

    def F(self, b, x):
        """Gradient (first derivative)"""
        alpha.vector()[:] = x
        assemble(E_alpha, b)

    def J(self, A, x):
        """Hessian (second derivative)"""
        alpha.vector()[:] = x
        assemble(E_alpha_alpha, A)
        
#PETScOptions.set("help")
solver_alpha_tao = PETScTAOSolver()
#solver_alpha_tao.parameters.update({"monitor_convergence": True})
PETScOptions.set("tao_type","tron")
#PETScOptions.set("pc_type","ilu")
#PETScOptions.set("ksp_type","nash")
PETScOptions.set("tao_monitor")
lb = interpolate(Constant("0."), V_alpha) # lower bound, initialize to 0
ub = interpolate(Constant("1."), V_alpha) # upper bound, set to 1
for bc in bc_alpha:
    bc.apply(lb.vector())
    bc.apply(ub.vector())
    
    
    
def alternate_minimization(u, alpha, tol=1.e-2, maxiter=100, alpha_0 = interpolate(Constant("0.0"), V_alpha)):
    # initialization
    iter = 1; err_alpha = 1
    alpha_error = Function(V_alpha)
    # iteration loop
    i_t = 1
    
    while err_alpha>tol and iter<maxiter:
        plt.figure(i_t)
        plt.colorbar(plot(alpha, range_min=0., range_max=1., title = "Damage at loading %.4f"%(t*load0)))
        plt.savefig(directory+'/alpha' + str(i_t) + '.png')
        plt.close(i_t)
        i_t += 1   
        # solve elastic problem
        solver_u.solve()
        # solve damage problem
        #solver_alpha.solve()
        solver_alpha_tao.solve(DamageProblem(), alpha.vector(), lb.vector(), ub.vector())# test error
        alpha_error.vector()[:] = alpha.vector() - alpha_0.vector()
        err_alpha = norm(alpha_error.vector(), "linf")
        # monitor the results
        if MPI.comm_world.rank == 0:
            print("Iteration:  %2d, Error: %2.8g, alpha_max: %.8g" %(iter, err_alpha, alpha.vector().max()))
        # update iteration
        alpha_0.assign(alpha)
        iter=iter+1
    return (err_alpha, iter)



##########################################################




    

def postprocessing():
    # plt.figure(i_t)
    # plt.colorbar(plot(alpha, range_min=0., range_max=1., title = "Damage at loading %.4f"%(t*load0)))
    # plt.savefig(savedir+'/alpha' + str(t) + '.png')
    # # Save number of iterations for the time step
    iterations[i_t] = np.array([t,i_t])
    # Calculate the energies
    elastic_energy_value = assemble(Pi)
    surface_energy_value = assemble(dissipated_energy)

    if MPI.comm_world.rank == 0:
        print("\nEnd of timestep %d with load multiplier %g"%(i_t, t))
        print("\nElastic and surface energies: (%g,%g)"%(elastic_energy_value,surface_energy_value))
        print("-----------------------------------------")
    energies[i_t] = np.array([t,elastic_energy_value,surface_energy_value,elastic_energy_value+surface_energy_value])
    # Calculate the axial force resultant
    forces[i_t] = np.array([t,Piola1st[0,0]((distance_between_voids+R_void,0)), Cauchy[0,0]((distance_between_voids+R_void,0))])
    # Dump solution to file
    # p = q.sub(1, True)
    # p.rename("pressure", "pressure")
    # output_file.write(p,t)
    u = q.sub(0, True)
    u.rename("displacement", "displacement")
    output_file.write(u,t)

    # forces[i_t] = np.array([t,Piola1st[0,0]((distance_between_voids+R_void,0)), p((distance_between_voids+R_void,0))])

    
    strain_energy = project(psi, V_alpha)
    strain_energy.rename("strain_energy","strain_energy")
    output_file.write(strain_energy, t)
    alpha.rename("damage", "damage")
    output_file.write(alpha,t)
    
    displacements[i_t] = np.array([t,u[0]((distance_between_voids+R_void,0)) , u[0]((1,0))])
    # Save some global quantities as a function of the time
    np.savetxt(savedir+'/energies.txt', energies)
    np.savetxt(savedir+'/forces.txt', forces)
    np.savetxt(savedir+'/iterations.txt', iterations)
    np.savetxt(savedir+'/displacements.txt', displacements)
    
    
################################################


lambdab = (R+load_multipliers)/R 

lambdaa = ((lambdab**2-1)*(R/R_void)**2 + 1)**0.5

#  loading and initialization of vectors to store time datas
energies = np.zeros((len(load_multipliers),4))
iterations = np.zeros((len(load_multipliers),2))

forces = np.zeros((len(load_multipliers),3))
displacements = np.zeros((len(load_multipliers),3))

lb.interpolate(Constant(0.))
for (i_t, t) in enumerate(load_multipliers):

    xmove.t = 1*t
    ymove.t = 1*t


    # solve alternate minimization
    solver_u.solve()
    # u = q.sub(0, True)
    # u.rename("displacement", "displacement")
    # output_file.write(u,t)

    # alternate_minimization(u,alpha,maxiter=30)
    
    u,p = split(q)

    DeformationGradient = grad(u) + I

    # Right Cauchy-Green tensor
    RCGT = DeformationGradient.T * DeformationGradient
    LCGT = DeformationGradient * DeformationGradient.T

    CauchyST = - p * I + mu * LCGT
    Piola1st = ufl.variable(diff(a(alpha)* psi  , F))

    Piola1st = project(Piola1st, V_tensor)
    Piola1st.rename("StressTensor", "StressTensor")
    output_file.write(Piola1st,t)
    print("External displacement:", u[0]((1,0)))

    Cauchy = project(CauchyST, V_tensor)
    print(Cauchy[0,0]((1,0)))

    Cauchy.rename("CauchyStressTensor", "CauchyStressTensor")
    output_file.write(Cauchy,t)
    
    # updating the lower bound to account for the irreversibility
    lb.vector()[:] = alpha.vector()
    postprocessing()


def plot_stretches():
    plt.figure()
    plt.plot(lambdab, lambdaa)
    plt.plot((displacements[:,2] + R)/R, (displacements[:,1]+R_void)/R_void, '*')
    plt.xlabel(r'\lambda_b')
    plt.ylabel(r'\lambda_a')
    plt.savefig(savedir+'/Stretches.png')
    # plt.show()

# plot_stretches()