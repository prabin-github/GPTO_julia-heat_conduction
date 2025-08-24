function init_geometry(FE,OPT,GEOM)
    #
    # Initialize GEOM structure with initial design
    #

    if !GEOM[:initial_design][:restart]
        # Run the initial design file (like exec in Python)
        include(GEOM[:initial_design][:path])

    # To use non contiguous numbers in the point_mat, we need to grab the
    # points whose ID matches the number specified by bar_mat. We achieve 
    # this via a map (sparse vector) between point_mat_rows and pt_IDs st
    # point_mat_row(point_ID) = row # of point_mat for point_ID
    pt_IDs = GEOM[:initial_design][:point_matrix][:, 1]  # Julia is 1-based

    GEOM[:point_mat_row] = sparse(pt_IDs)
    else
        include(GEOM[:initial_design][:path])
        GEOM[:initial_design][:point_matrix] = GEOM[:current_design][:point_matrix]
        GEOM[:initial_design][:bar_matrix] = GEOM[:current_design][:bar_matrix]
    end

    GEOM[:n_point] = size(GEOM[:initial_design][:point_matrix], 1)
    GEOM[:n_bar] = size(GEOM[:initial_design][:bar_matrix], 1)
end


function compute_bar_elem_distance(FE,OPT,GEOM)
    #
    # This function computes an array dist of dimensions n_bar x n_elem with 
    # the distance from the centroid of each element to each bar's medial axis.
    # 
    # Ddist_Dbar_ends is a 3-dimensional array of dimensions n_bar_dofs x n_bar
    # x n_elem that contains the sensitivities of the signed distances in dist
    # with respect to each of the n__bar_dofs coordinates of the medial axis 
    # end points.

    # global FE, GEOM, OPT

    ## set parameters
    tol = 1e-12 # tolerance on the length of a bar

    n_elem = FE[:n_elem]
    dim = FE[:dim]
    n_bar = GEOM[:n_bar]
    n_bar_dofs = 2 * dim

    # The following code is vectorized over the bars and the elements. We
    # have to be consistent with the order of indices to perform element-
    # wise array operations. The order of indices is (dim,bar,element)
    points = transpose(GEOM[:current_design][:point_matrix][:, 2:end])

    x_1b = points[Int.(OPT[:bar_dv][collect(1:dim), :])]           # (i,b) 
    x_2b = points[Int.(OPT[:bar_dv][collect(dim .+ 1:2*dim), :])]  # (i,b) 
    x_e = reshape(FE[:centroids], (size(FE[:centroids], 1), 1, size(FE[:centroids], 2)))  # (i,1,e) 

    a_b = x_2b - x_1b                   # Numerator of Eq. (11)
    l_b = sqrt.(sum(a_b .^ 2, dims=1))  # length of the bars, Eq. (10)
    l_b[vec(l_b .< tol)] .= 1           # To avoid division by zero
    a_b = a_b ./ l_b                    # normalize the bar direction to unit vector, Eq. (11)

    x_e_1b = x_e .- x_1b                 # (i,b,e) 
    x_e_2b = x_e .- x_2b                 # (i,b,e) 
    norm_x_e_1b = sqrt.(sum(x_e_1b .^ 2, dims=1)) # (1,b,e)
    norm_x_e_2b = sqrt.(sum(x_e_2b .^ 2, dims=1)) # (1,b,e) 

    l_be = sum(x_e_1b .* a_b, dims=1)        # (1,b,e), Eq. (12)
    vec_r_be = x_e_1b - l_be .* a_b          # (i,b,e)
    r_be = sqrt.(sum(vec_r_be .^ 2, dims=1)) # (1,b,e), Eq. (13)

    branch1 = l_be .<= 0.0            # (1,b,e)
    branch2 = l_be .> l_b             # (1,b,e)
    branch3 = .~(branch1 .| branch2)  # (1,b,e)

    # Compute the distances, Eq. (14)
    dist_tmp = branch1 .* norm_x_e_1b .+
               branch2 .* norm_x_e_2b .+
               branch3 .* r_be               # (1,b,e)
    dist = permutedims(dist_tmp, (2, 3, 1))  # (b,e)

    ## compute sensitivities
    Dd_be_Dx_1b = zeros(FE[:dim], n_bar, n_elem)
    Dd_be_Dx_2b = zeros(FE[:dim], n_bar, n_elem)

    d_inv = dist_tmp .^ -1    # This can render a division by zero (if point 
    d_inv[isinf.(d_inv)] .= 0 # lies on medial axis), and so we now fix it
    l_be_over_l_b = l_be ./ l_b

    ## The sensitivities below are obtained from Eq. (30)
    ## sensitivity to x_1b
    if sum(branch1) > 0
        branch1_1 = repeat(branch1, 2, 1, 1)
        d_inv_1 = repeat(d_inv, 2, 1, 1)
        Dd_be_Dx_1b[branch1_1] = -1 * x_e_1b[branch1_1] .* d_inv_1[branch1_1]
    end
    # Dd_bd_Dx_1b(:,branch2) = 0;
    if sum(branch3) > 0
        branch3_1 = repeat(branch3, 2, 1, 1)
        d_inv_1 = repeat(d_inv, 2, 1, 1)
        l_be_over_l_b_1 = repeat(l_be_over_l_b, 2, 1, 1)
        Dd_be_Dx_1b[branch3_1] = -vec_r_be[branch3_1] .* d_inv_1[branch3_1] .* (1 .- l_be_over_l_b_1[branch3_1])
    end

    ## sensitivity to x_2b
    # Dd_bd_Dx_2b(:,branch1) = 0;
    if sum(branch2) != 0
        branch2_1 = repeat(branch2, 2, 1, 1)
        d_inv_1 = repeat(d_inv, 2, 1, 1)
        Dd_be_Dx_2b[branch2_1] = -x_e_2b[branch2_1] .* d_inv_1[branch2_1]
    end
    if sum(branch3) != 0
        branch3_1 = repeat(branch3, 2, 1, 1)
        d_inv_1 = repeat(d_inv, 2, 1, 1)
        l_be_over_l_b_1 = repeat(l_be_over_l_b, 2, 1, 1)
        Dd_be_Dx_2b[branch3_1] =
            -vec_r_be[branch3_1] .* d_inv_1[branch3_1] .* l_be_over_l_b_1[branch3_1]
    end

    ## assemble the sensitivities to the bar design parameters (scaled)
    Ddist_Dbar_ends = zeros(n_bar_dofs, n_bar, n_elem)
    Ddist_Dbar_ends[collect(1:dim), :, :]= Dd_be_Dx_1b .* OPT[:scaling][:point_scale]
    Ddist_Dbar_ends[dim .+ collect(1:dim), :, :] = Dd_be_Dx_2b .* OPT[:scaling][:point_scale]

    Ddist_Dbar_ends = permutedims(Ddist_Dbar_ends, (2, 3, 1))

    return dist, Ddist_Dbar_ends
end


function penalize(args...)
    # [P, dPdx] = penalize(x, p, penal_scheme)
    #     penalize(x) assumes x \in [0,1] and decreases the intermediate values
    #
    #	  For a single input, the interpolation is SIMP with p = 3
    #
    #	  The optional second argument is the parameter value p.
    #
    #     The optional third argument is a string that indicates the way the 
    #	  interpolation is defined, possible values are:
    #       'SIMP'      : default 
    # 	  	'RAMP'      : 
    #

    # consider input
    n_inputs = length(args)
    x = args[1]
    if n_inputs == 1
        # set the definition to be used by default.
        p = 3
        penal_scheme = "SIMP"
    elseif n_inputs == 2
        p = args[2]
        penal_scheme = "SIMP"
    elseif n_inputs == 3
        p = args[2]
        penal_scheme = args[3]
    end

    # consider output
    ### not implemented

    # define
    if penal_scheme == "SIMP"
        P    = x .^ p
        dPdx = p .* x .^ (p - 1)
    elseif penal_scheme == "RAMP"
        P    = x ./ (1 .+ p .* (1 .- x))
        dPdx = (1 .+ p) ./ (1 .+ p .* (1 .- x)).^2
    else
        println("Unidentified parameters")
        return nothing, nothing
    end

    return P, dPdx
end


function project_element_densities(FE,OPT,GEOM)
    # This def computes the combined unpenalized densities (used to
    # compute the volume) and penalized densities (used to compute the ersatz
    # material for the analysis) and saves them in the global variables
    # FE['elem_dens'] and FE['penalized_elem_dens'].  
    #
    # It also computes the derivatives of the unpenalized and penalized
    # densities with respect to the design parameters, and saves them in the
    # global variables FE['Delem_dens_Ddv'] and FE['Dpenalized_elem_dens_Ddv']. 
    #

    ##  Distances from the element centroids to the medial segment of each bar
    d_be, Dd_be_Dbar_ends = compute_bar_elem_distance(FE,OPT,GEOM)

    ## Bar-element projected densities
    r_b = GEOM[:current_design][:bar_matrix][:, end] # bar radii
    r_e = OPT[:parameters][:elem_r] # sample window radius
    
    # X_be is \phi_b/r in Eq. (2).  Note that the numerator corresponds to
    # the signed distance of Eq. (8).
    X_be = (r_b .- d_be[:, :, 1]) ./ transpose(r_e)

    # Projected density 
    if FE[:dim] == 2  # 2D
        rho_be = ((π .- acos.(Complex.(X_be)) .+ X_be .* sqrt.(-1 * Complex.(X_be) .^ 2 .+ 1.0)) ./ π) .* (abs.(Complex.(X_be)) .< 1) .+ 1 * (X_be .>= 1)
        Drho_be_Dx_be = ((sqrt.(-(Complex.(X_be) .^ 2) .+ 1.0) .* 2.0) ./ π) .* (abs.(Complex.(X_be)) .< 1) # Eq. (28)
    elseif FE[:dim] == 3
        rho_be = (((X_be .+ 1.0) .^ 2) .* (X_be .- 2.0) .* (-1.0 ./ 4.0)) .* (abs.(X_be) .< 1) .+ 1 * (X_be .>= 1) # Eqs. (2) and (3)
        Drho_be_Dx_be = ((X_be .^ 2) .* (-3.0 ./ 4.0) .+ 3.0 ./ 4.0) .* (abs.(X_be) .< 1)# Eq. (28)
    end
    rho_be = real(rho_be)
    Drho_be_Dx_be = real(Drho_be_Dx_be)

    # Sensitivities of raw projected densities, Eqs. (27) and (29)
    Drho_be_Dbar_ends = Drho_be_Dx_be .* (-transpose(r_e) .^ (-1)) .* Dd_be_Dbar_ends
    Drho_be_Dbar_radii = Drho_be_Dx_be .* (transpose(r_e) .^ (-1)) .* OPT[:scaling][:radius_scale]

    ## Combined densities
    # Get size variables
    alpha_b = GEOM[:current_design][:bar_matrix][:, end-1] # bar size

    # Without penalization:
    # ====================
    # X_be here is \hat{\rho}_b in Eq. (4) with the value of q such that
    # there is no penalization (e.g., q = 1 in SIMP).
    X_be = rho_be .* alpha_b

    # Sensitivities of unpenalized effective densities, Eq. (26) with
    # ?\partial \mu / \partial (\alpha_b \rho_{be})=1
    DX_be_Dbar_ends = Drho_be_Dbar_ends .* alpha_b
    DX_be_Dbar_size = copy(rho_be)
    DX_be_Dbar_radii = Drho_be_Dbar_radii .* alpha_b

    # Combined density of Eq. (5).
    rho_e, Drho_e_DX_be = smooth_max(X_be,
        OPT[:parameters][:smooth_max_param],
        OPT[:parameters][:smooth_max_scheme],
        FE[:material][:rho_min])

    # Sensitivities of combined densities, Eq. (25)
    Drho_e_Dbar_ends = Drho_e_DX_be .* DX_be_Dbar_ends
    Drho_e_Dbar_size = Drho_e_DX_be .* DX_be_Dbar_size
    Drho_e_Dbar_radii = Drho_e_DX_be .* DX_be_Dbar_radii

    # Stack together sensitivities with respect to different design
    # variables into a single vector per element
    Drho_e_Ddv = zeros(FE[:n_elem], OPT[:n_dv])
    for b = 1:GEOM[:n_bar]
        Drho_e_Ddv[:, Int.(OPT[:bar_dv][:, b])] = Drho_e_Ddv[:, Int.(OPT[:bar_dv][:, b])] .+
                                                hcat(reshape(Drho_e_Dbar_ends[b, :, :], (FE[:n_elem], 2 * FE[:dim])),
            reshape(Drho_e_Dbar_size[b, :], (FE[:n_elem], 1)),
            reshape(Drho_e_Dbar_radii[b, :], (FE[:n_elem], 1)))
    end

    # With penalization:   
    # =================
    # In this case X_be *is* penalized (Eq. (4)).
    penal_X_be, Dpenal_X_be_DX_be = penalize(X_be,
        OPT[:parameters][:penalization_param],
        OPT[:parameters][:penalization_scheme])

    # Sensitivities of effective (penalized) densities, Eq. (26)
    Dpenal_X_be_Dbar_ends = Dpenal_X_be_DX_be .* DX_be_Dbar_ends
    Dpenal_X_be_Dbar_size = Dpenal_X_be_DX_be .* DX_be_Dbar_size
    Dpenal_X_be_Dbar_radii = Dpenal_X_be_DX_be .* DX_be_Dbar_radii

    # Combined density of Eq. (5).
    penal_rho_e, Dpenal_rho_e_Dpenal_X_be = smooth_max(penal_X_be,
            OPT[:parameters][:smooth_max_param],
            OPT[:parameters][:smooth_max_scheme],
            FE[:material][:rho_min])

    # Sensitivities of combined densities, Eq. (25)
    Dpenal_rho_e_Dbar_ends = Dpenal_rho_e_Dpenal_X_be .* Dpenal_X_be_Dbar_ends
    Dpenal_rho_e_Dbar_size = Dpenal_rho_e_Dpenal_X_be .* Dpenal_X_be_Dbar_size
    Dpenal_rho_e_Dbar_radii = Dpenal_rho_e_Dpenal_X_be .* Dpenal_X_be_Dbar_radii

    # Sensitivities of projected density
    Dpenal_rho_e_Ddv = zeros(FE[:n_elem], OPT[:n_dv])
    
    # Stack together sensitivities with respect to different design
    # variables into a single vector per element
    for b = 1:GEOM[:n_bar]
        Dpenal_rho_e_Ddv[:, Int.(OPT[:bar_dv][:, b])] = Dpenal_rho_e_Ddv[:, Int.(OPT[:bar_dv][:, b])] +
                                                      hcat(reshape(Dpenal_rho_e_Dbar_ends[b, :, :], (FE[:n_elem], 2 * FE[:dim])),
            reshape(Dpenal_rho_e_Dbar_size[b, :], (FE[:n_elem], 1)),
            reshape(Dpenal_rho_e_Dbar_radii[b, :], (FE[:n_elem], 1)))
    end

    ## Write the element densities and their sensitivities to OPT 
    OPT[:elem_dens] = copy(rho_e[:])
    OPT[:Delem_dens_Ddv] = copy(Drho_e_Ddv)
    OPT[:penalized_elem_dens] = copy(penal_rho_e[:])
    OPT[:Dpenalized_elem_dens_Ddv] = copy(Dpenal_rho_e_Ddv)
end


function smooth_max(x,p,form_def,x_min)
    #
    # This def computes a smooth approximation of the maximum of x.  The
    # type of smooth approximation (listed below) is given by the argument
    # form_def, and the corresponding approximation parameter is given by p.
    # x_min is a lower bound to the smooth approximation for the modified
    # p-norm and modified p-mean approximations.
    #
    #
    #     The optional third argument is a string that indicates the way the 
    #	  approximation is defined, possible values are:
    # 		'mod_p-norm'   : overestimate using modified p-norm (supports x=0)
    # 		'mod_p-mean'   : underestimate using modified p-norm (supports x=0)
    #		'KS'           : Kreisselmeier-Steinhauser, overestimate
    #		'KS_under'     : Kreisselmeier-Steinhauser, underestimate
    #

    if form_def== "mod_p-norm"
        # Eq. (6)
        # in this case, we assume x >= 0 
        S = ( x_min.^p .+ (1 .- x_min.^p)*sum(x.^p,dims=1) ).^(1/p)
        dSdx = (1-x_min^p)*(x./S).^(p-1)
    elseif form_def == "mod_p-mean"
        # in this case, we assume x >= 0 
        N = size(x, 1)
        S = (x_min^p + (1 - x_min^p) * sum(x.^p, dims=1) / N) .^ (1/p)
        dSdx = (1 - x_min^p) * (1/N) * (x ./ S) .^ (p - 1)
    elseif form_def == "KS"
        epx = exp.(p .* x)
        sum_epx = sum(epx, dims=1)
        S = x_min .+ (1 - x_min) .* log.(sum_epx) ./ p
        dSdx = (1 - x_min) .* epx ./ sum_epx
    elseif form_def == "KS_under"
        N = size(x, 1)
        epx = exp.(p .* x)
        sum_epx = sum(epx, dims=1)
        S = x_min .+ (1 - x_min) .* log.(sum_epx ./ N) ./ p
        dSdx = (1 - x_min) .* epx ./ sum_epx
    ### softmax added ###
    elseif form_def == "softmax"
        epx     = exp.(p .* (x .- maximum(x)))
        sum_epx = sum(epx, dims=1)
        S       = sum(x .* epx, dims=1) ./ sum_epx
        dSdx    = (epx ./ sum_epx) .* (1 .+ p .* (x .- S))
    ### modified softmax added ###
    elseif form_def == "mod_softmax"
        epx     = exp.(p .* (x .- maximum(x)))
        sum_epx = sum(epx, dims=1)
        A       = sum(x .* epx, dims=1) ./ sum_epx
        S       = x_min^p .+ (1 - x_min^p) * A
        dSdx    = (1 - x_min^p) .* (epx ./ sum_epx) .* (1 .+ p .* (x .- A))
    else
        println("\nsmooth_max received invalid form_def.\n")
        return nothing, nothing
    end

    return S, dSdx
end


function update_dv_from_geom(FE,OPT,GEOM)
    #
    # This def updates the values of the design variables (which will be
    # scaled if OPT.options['dv']_scaling is true) based on the unscaled bar 
    # geometric parameters. It does the opposite from the def 
    # update_geom_from_dv.
    #

    # global GEOM, OPT 

    # Fill in design variable vector based on the initial design
    # Eq. (32)
    OPT[:dv][OPT[:point_dv]] = (transpose(GEOM[:initial_design][:point_matrix][:, 2:end]) .- 
                OPT[:scaling][:point_min]) ./ OPT[:scaling][:point_scale]
    OPT[:dv][OPT[:size_dv]] = GEOM[:initial_design][:bar_matrix][:, end-1]
    OPT[:dv][Int.(OPT[:radius_dv])] = (GEOM[:initial_design][:bar_matrix][:, end] .- 
                OPT[:scaling][:radius_min]) ./ OPT[:scaling][:radius_scale]
end


function update_geom_from_dv(FE,OPT,GEOM)
    # This def updates the values of the unscaled bar geometric parameters
    # from the values of the design variableds (which will be scaled if
    # OPT.options['dv']_scaling is true). It does the
    # opposite from the def update_dv_from_geom.
    #
    # global GEOM , OPT , FE

    # Eq. (32)
    GEOM[:current_design][:point_matrix][:, 2:end] = transpose(OPT[:scaling][:point_scale] .* 
        reshape(OPT[:dv][OPT[:point_dv]], (FE[:dim], GEOM[:n_point])) .+ OPT[:scaling][:point_min])
    
    GEOM[:current_design][:bar_matrix][:, end-1] = OPT[:dv][OPT[:size_dv]]
    
    GEOM[:current_design][:bar_matrix][:, end] = OPT[:dv][Int.(OPT[:radius_dv])] .* 
        OPT[:scaling][:radius_scale] .+ OPT[:scaling][:radius_min]
end