function FE_solve(FE,OPT,GEOM)
    # This function solves the system of linear equations arising from the
    # finite element discretization of Eq. (17).  It stores the displacement 
    # and the reaction forces in FE['U'] and FE['P'].

    f(x) = 1.0     # source term (uniform heat generation)
    g(x) = 0.0     # zero temperature for Dirichlet BC
    h(x) = 0.0     # adiabatic boundary condition for Neumann BC

    # Constitutive equation (Fourier's Law)
    ϕ(ρ) = FE[:material][:Kt_min] + ρ*(FE[:material][:Kt]-FE[:material][:Kt_min])
    q(ρ,∇u) = ϕ(ρ) * ∇u
    dq(ρ,∇u) = ϕ'(ρ) * ∇u

    # === FE analysis ===
    # Weak bilinear form
    a(u,v) = ∫( ∇(v) ⊙ ( q∘(OPT[:penalized_elem_dens],∇(u)) ) ) * GEOM[:param].dΩ
    b(v) = ∫( v*f ) * GEOM[:param].dΩ + ∫( v*h ) * GEOM[:param].dΓ
    op = AffineFEOperator(a, b, GEOM[:param].U, GEOM[:param].V0)
    FE[:uh] = solve(op)
end


function init_FE(FE,OPT,GEOM)
    # Initialize the FE structure
    # global FE

    if FE[:mesh_input][:type] == "generate"
        generate_mesh(FE,OPT,GEOM)
    # elseif FE[:mesh_input][:type] == "read-home-made"
    #     load(FE[:mesh_input][:mesh_filename])
    elseif FE[:mesh_input][:type] == "read-gmsh"
        read_gmsh(FE,OPT,GEOM)
    else
        println("Unidentified mesh type")
    end

    # Setup boundary conditions
    include(FE[:mesh_input][:bcs_file])
end
