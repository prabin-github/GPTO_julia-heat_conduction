## Input file 
#
# *** THIS SCRIPT HAS TO BE CUSTOMIZED BY THE USER ***
#
# This script sets up the displacement boundary conditions and the forces
# for the analysis. 
#
# Important note: you must make sure you do not simultaneously impose 
# displacement boundary conditions and forces on the same degree of
# freedom.

# ** Do not modify this line **

nelx = Int(FE[:mesh_input][:elements_per_side][1])
nely = Int(FE[:mesh_input][:elements_per_side][2])
## ============================        


# The model can be written to a vtk file and opened with ParaView to 
# determine the tags of different entities, which we need to impose BCs
labels = get_face_labeling(GEOM[:model])
################## DOFs for Dirichlet BC (Zero temperature) ##################
entity = num_entities(labels)+1
mid = fld(nelx, 2)+1
span = round(Int, nelx/20)
fixeddofs_nodes = collect(mid-span:mid+span)      # L/10 mid-nodes at bottom edge for Dirichlet BC
labels.d_to_dface_to_entity[1][fixeddofs_nodes] .= entity
add_tag!(labels,"fixeddofs_nodes",[entity])
##############################################################################
# writevtk(model,"heat_conduction")


# Create FE spaces
order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V0 = TestFESpace(GEOM[:model],reffe;
  conformity=:H1,
  dirichlet_tags="fixeddofs_nodes")

f(x) = 1.0     # source term (uniform heat generation)
g(x) = 0.0     # zero temperature for Dirichlet BC
h(x) = 0.0     # adiabatic boundary condition for Neumann BC
U = TrialFESpace(V0, g)

# Now extract DOF indices
free_ids = get_free_dof_ids(U)          # indices of unknown DOFs
fixed_ids  = get_dirichlet_dof_ids(U)   # indices of fixed (Dirichlet) DOFs

# Extract the Dirichlet‐DOF values from the FESpace
fixeddofs_values = get_dirichlet_dof_values(U)

# Integration mesh
degree = 2*order
Ω = Triangulation(GEOM[:model])
dΩ = Measure(Ω,degree)

# Boundary triangulation for Neumann BCs
neumanntags = ["boundary"]
Γ = BoundaryTriangulation(GEOM[:model],tags=neumanntags)
dΓ = Measure(Γ,degree)

GEOM[:param] = (; V0, U, Ω, dΩ, Γ, dΓ)

FE[:BC] = Dict{Symbol, Any}()
FE[:BC][:n_pre_disp_dofs] = length(fixeddofs_nodes)
FE[:BC][:disp_node] = fixeddofs_nodes
FE[:BC][:disp_value] = fixeddofs_values
