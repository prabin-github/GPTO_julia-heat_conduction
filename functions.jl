function compute_compliance(FE,OPT,GEOM)
    #
    # This function computes the mean compliance and its sensitivities
    # based on the last finite element analysis
    # global FE, OPT

    # Constitutive equation (Fourier's Law)
    ϕ(ρ) = FE[:material][:Kt_min] + ρ*(FE[:material][:Kt]-FE[:material][:Kt_min])
    q(ρ,∇u) = ϕ(ρ) * ∇u
    dq(ρ,∇u) = ϕ'(ρ) * ∇u

    # === Objective function and sensitivity analysis ===
    c = sum( ∫( ∇(FE[:uh]) ⊙ ( q∘(OPT[:penalized_elem_dens],∇(FE[:uh] ) ) ) ) * GEOM[:param].dΩ )
    dc = ∫( -∇(FE[:uh]) ⊙ ( dq∘(OPT[:penalized_elem_dens],∇(FE[:uh])) ) ) * GEOM[:param].dΩ
    dc = vec(reshape(get_contribution(dc, GEOM[:param].Ω), (1, FE[:n_elem])) * OPT[:Dpenalized_elem_dens_Ddv])

    OPT[:compliance] = c
    OPT[:grad_compliance] = dc

    return c, dc
end


function compute_volume_fraction(FE,OPT,GEOM)
    #
    # This function computes the volume fraction and its sensitivities
    # based on the last geometry projection
    #
    # global FE, OPT

    # compute the volume fraction
    v_e = FE[:elem_vol] # element volume
    V = sum(v_e)         # full volume
    v = dot(v_e, OPT[:elem_dens]) # projected volume
    volfrac = v / V # Eq. (16)

    # compute the design sensitivity
    Dvolfrac_Ddv = (transpose(v_e[:]) * OPT[:Delem_dens_Ddv]) / V   # Eq. (31)
    grad_vofrac = transpose(Dvolfrac_Ddv)

    # output
    OPT[:volume_fraction] = volfrac
    OPT[:grad_volume_fraction] = grad_vofrac

    return volfrac, grad_vofrac
end


function evaluate_relevant_functions(FE,OPT,GEOM)
    # Evaluate_relevant_functions() looks at OPT['functions'] and evaluates the
    # relevant functions for this problem based on the current OPT['dv']
    # global OPT

    OPT[:functions][:n_func] = length(OPT[:functions][:f])

    for i in 1:OPT[:functions][:n_func]
        funcname = OPT[:functions][:f][i][:function]
        value, grad = eval(Meta.parse("$(funcname)(FE, OPT, GEOM)"))
        OPT[:functions][:f][i][:value] = value
        OPT[:functions][:f][i][:grad] = grad
    end
end


function nonlcon(FE,OPT,GEOM,dv)
    # global FE, OPT, GEOM

    OPT[:dv_old] = copy(OPT[:dv])
    OPT[:dv] = copy(dv)

    # Update geometry and analysis only if design variables changed
    if any(OPT[:dv] .!= OPT[:dv_old])
        update_geom_from_dv(FE,OPT,GEOM)
        perform_analysis(FE,OPT,GEOM)
    end

    n_con = OPT[:functions][:n_func]-1  # number of constraints
    g = zeros(n_con)

    for i in 1:n_con
        g[i] = OPT[:functions][:f][i+1][:value]
    end
    g .-= OPT[:functions][:constraint_limit]

    return g
end


function nonlcongrad(FE,OPT,GEOM,dv)
    # global FE, OPT, GEOM

    n_con = OPT[:functions][:n_func]-1   # number of constraints (exclude the objective)
    gradg = zeros(Float64, OPT[:n_dv], n_con)

    for i in 1:n_con
        gradg[:, i] = OPT[:functions][:f][i+1][:grad]
    end
    
    return vec(gradg)
end


function obj(FE,OPT,GEOM,dv)
    # global  FE, OPT, GEOM

    OPT[:dv_old] = copy(OPT[:dv])   # save the previous design
    OPT[:dv] = copy(dv)              # update the design

    # If different, update or perform the analysis
    if any(OPT[:dv] .!= OPT[:dv_old])
        update_geom_from_dv(FE,OPT,GEOM)
        perform_analysis(FE,OPT,GEOM)
    end

    f = copy(OPT[:functions][:f][1][:value])
    g = copy(OPT[:functions][:f][1][:grad])
    
    return f, g
end


function perform_analysis(FE,OPT,GEOM)
    # Perform the geometry projection, solve the finite
    # element problem for the displacements and reaction forces, and then
    # evaluate the relevant functions.
    project_element_densities(FE,OPT,GEOM)
    FE_solve(FE,OPT,GEOM)    
    evaluate_relevant_functions(FE,OPT,GEOM)
end
