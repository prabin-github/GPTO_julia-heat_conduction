function init_optimization(FE,OPT,GEOM)
    # Initialize functions to compute
    # Concatenate list of functions to be computed
    f_list = Dict{Int, String}()
    f_list[1] = OPT[:functions][:objective]
    f_list[2] = OPT[:functions][:constraints]

    # here we list all the functions that are available to compute as f{i}
    f = Dict{Int, Dict{Symbol, String}}()

    f[1] = Dict(
        :name  => "compliance",
        :function  => "compute_compliance")

    f[2] = Dict(
        :name  => "volume fraction",
        :function  => "compute_volume_fraction")

    # compare all functions available with the ones specified in inputs.m
    n = length(f)
    m = length(f_list)

    OPT[:functions][:f] = Dict{Int, Dict{Symbol, Any}}()
    for j in 1:m
        for i in 1:n
            if f[i][:name] == f_list[j]
                OPT[:functions][:f][j] = f[i]
            end
        end
    end

    OPT[:functions][:n_func] = length(OPT[:functions][:f])

    ## initialize sample window size
    if !haskey(OPT[:parameters], :elem_r)
        # compute sampling radius
        # The radius corresponds to the circle (or sphere) that circumscribes a
        # square (or cube) that has the edge length of elem_size.
        OPT[:parameters][:elem_r] = sqrt(FE[:dim])/2 * FE[:elem_vol].^(1/FE[:dim])
    end

    ##
    # Initilize the design variable and its indexing schemes

    # we are designing the points, the size variables, and the radii of the
    # bars:
    OPT[:n_dv] = FE[:dim] * GEOM[:n_point] + 2 * GEOM[:n_bar]
    OPT[:dv] = zeros(OPT[:n_dv])

    OPT[:point_dv] = collect(1:FE[:dim]*GEOM[:n_point]) # such that dv(point_dv) = point
    OPT[:size_dv] = OPT[:point_dv][end] .+ collect(1:GEOM[:n_bar])
    OPT[:radius_dv] = OPT[:size_dv][end] .+ collect(1:GEOM[:n_bar])

    OPT[:scaling] = Dict{Symbol, Any}()

    if OPT[:options][:dv_scaling] == 1
        # Compute variable limits for Eq. (32)
        OPT[:scaling][:point_scale] = FE[:coord_max] - FE[:coord_min]
        OPT[:scaling][:point_min] = FE[:coord_min]
        # Consider possibility that max_bar_radius and min_bar_radius are
        # the same (when bars are of fixed radius)
        delta_radius = GEOM[:max_bar_radius] - GEOM[:min_bar_radius]
        if delta_radius < 1e-12
            OPT[:scaling][:radius_scale] = 1
        else
            OPT[:scaling][:radius_scale] = delta_radius
        end
        OPT[:scaling][:radius_min] = GEOM[:min_bar_radius]
    else
        OPT[:scaling][:point_scale] = 1.0
        OPT[:scaling][:point_min] = 0.0
        OPT[:scaling][:radius_scale] = 1.0
        OPT[:scaling][:radius_min] = 0.0
    end

    # fill in design variable vector based on the initial design
    update_dv_from_geom(FE,OPT,GEOM)

    # set the current design to the initial design:
    GEOM[:current_design] = Dict{Symbol, Any}()
    GEOM[:current_design][:point_matrix] = GEOM[:initial_design][:point_matrix]
    GEOM[:current_design][:bar_matrix] = GEOM[:initial_design][:bar_matrix]

    # consider the bar design variables
    # Extract index of first and secont point of each bar
    x_1b_id = GEOM[:current_design][:bar_matrix][:, 2]
    x_2b_id = GEOM[:current_design][:bar_matrix][:, 3]

    # Extract index of first (second) point of each matrix
    pt1 = GEOM[:point_mat_row][x_1b_id]
    pt2 = GEOM[:point_mat_row][x_2b_id]

    pt_dv = reshape(OPT[:point_dv], (FE[:dim], GEOM[:n_point]))

    OPT[:bar_dv] = vcat(pt_dv[:, Int.(pt1)],
                        pt_dv[:, Int.(pt2)],
                        transpose(OPT[:size_dv]),
                        transpose(OPT[:radius_dv]))
end


function runmma(FE,OPT,GEOM,x0,obj,nonlcon)
    # 
    # Perform the optimization using MMA
    #
    # global OPT GEOM FE

    function plotfun(iter)
        if OPT[:options][:plot] == true
            plot_design(iter)

            if FE[:dim] == 2
                plot_density(iter)
                # plot_density_levelsets()
            end
        end
    end

    # Initialize history dictionary
    history = Dict{Symbol, Any}()

    # Design variables constraint
    # Initialize lower and upper bounds vectors
    if OPT[:options][:dv_scaling] == 1   # Eq. (33)
        lb_point = zeros(FE[:dim], 1)
        ub_point = ones(FE[:dim], 1)
        lb_radius = 0

        # Consider case when max_bar_radius and min_bar_radius are
        # the same (when bars are of fixed radius)
        if GEOM[:max_bar_radius] - GEOM[:min_bar_radius] < 1e-12
            ub_radius = 0
        else
            ub_radius = 1
        end
    else
        lb_point = FE[:coord_min]            # Eq. (18)
        ub_point = FE[:coord_max]            # Eq. (18)
        lb_radius = GEOM[:min_bar_radius]    # Eq. (19)
        ub_radius = GEOM[:max_bar_radius]    # Eq. (19)
    end

    lb_size = 0    # Eq. (20)
    ub_size = 1    # Eq. (20)

    lb_bar = [lb_point; lb_point; lb_size; lb_radius]
    ub_bar = [ub_point; ub_point; ub_size; ub_radius]

    lb = zeros(size(OPT[:dv]))
    ub = zeros(size(OPT[:dv]))

    lb[Int.(OPT[:bar_dv])] = repeat(lb_bar, 1, GEOM[:n_bar])
    ub[Int.(OPT[:bar_dv])] = repeat(ub_bar, 1, GEOM[:n_bar])

    ncons = OPT[:functions][:n_func]-1  # Number of optimization constraints
    ndv = copy(OPT[:n_dv]) # Number of design variables

    # Initialize vectors that store current and previous two design iterates
    x = copy(x0)
    xold1 = copy(x0)
    xold2 = copy(x0)

    # Initialize move limits
    ml_step = OPT[:options][:move_limit] * abs.(ub .- lb)  # Compute move limits once

    # Initialize lower and upper asymptotes
    low = copy(lb)
    upp = copy(ub)

    # These are the MMA constants (Svanberg, 1998 DACAMM Course)
    c = 1000*ones(ncons)
    d = ones(ncons)
    a0 = 1.0
    a = zeros(ncons)

    # Evaluate the initial design and print values to screen 
    iter = 0
    f0val, df0dx = obj(FE,OPT,GEOM,x)
    fval = nonlcon(FE,OPT,GEOM,x)[1]
    dfdx = nonlcongrad(FE,OPT,GEOM,x)

    dfdx = transpose(dfdx)      # m x n matrix to pass as dfdx into MMA

    # initialize history
    history[:fval] = Vector{Float64}()
    history[:fconsval] = Vector{Float64}()
    history[:x] = Vector{Vector{Float64}}()

    # Save history
    push!(history[:fval], f0val)
    push!(history[:fconsval], fval)
    push!(history[:x], x)
    
    #### Initialize stopping values
    kktnorm = 10*OPT[:options][:kkt_tol]
    dv_step_change = 10*OPT[:options][:step_tol]

    # Make sure the output folder exists, and if not, create it
    if !isdir(OPT[:options][:vtk_output_path])
        mkpath(OPT[:options][:vtk_output_path])
    end

    OPT[:options][:plot_output_path] = joinpath(OPT[:options][:vtk_output_path], "plot_design")
    
    if isdir(OPT[:options][:vtk_output_path])
        mkpath(OPT[:options][:plot_output_path])
    end

    # Produce output to screen
    filename = joinpath(OPT[:options][:vtk_output_path], OPT[:problem] * "_" * OPT[:functions][:objective] * "_" * 
                OPT[:functions][:constraints] * "_History.txt")

    fid = open(filename, "w")   # Open the file for writing (overwritten)
    println(fid, "Start time: ", Dates.format(now(), "yyyy-mm-dd E HH:MM:SS"))   # Write header
    println(fid, "It. ", iter, ", Obj = ", @sprintf("%.5e", f0val), ", ConsViol = ", @sprintf("%.5e", fval[1]), 
                    ", KKT-norm = ", @sprintf("%.5e", kktnorm), ", DV norm change = ", @sprintf("%.5e", dv_step_change))
    
    # Produce output to screen
    @printf("It. %d, Obj = %.5e, ConsViol = %.5e, KKT-norm = %.5e, DV norm change = %.5e\n",
            iter, f0val, fval[1], kktnorm, dv_step_change)

    # optional VTK output
    if OPT[:options][:write_to_vtk] == "all"
        writevtk(OPT[:options][:vtk_output_path], "dens", iter)
    end

    # Plot the initial design
    plotfun(iter)
    
    #
    # ******* MAIN MMA LOOP STARTS *******
    #
    while kktnorm > OPT[:options][:kkt_tol] && iter < OPT[:options][:max_iter] && 
    dv_step_change > OPT[:options][:step_tol]

        OPT[:iter] = iter
        iter = iter + 1

        ### Scaling of objective and its sensitivities ###
        f0val = f0val/history[:fval][1]
        df0dx = df0dx/history[:fval][1]
        ##################################################

        # Impose move limits by modifying lower and upper bounds passed to MMA, Eq. (33)
        mlb = max.(lb, x .- ml_step)
        mub = min.(ub, x .+ ml_step)

        #### Solve MMA subproblem for current design x
        xmma, ymma, zmma, lam, xsi, eta, mu, zet, s, low, upp =
            mmasub(ncons, ndv, iter, x, mlb, mub, xold1, xold2,
                f0val, df0dx, fval, dfdx, low, upp, a0, a, c, d)

        #### Updated design vectors of previous and current iterations
        xold2 .= copy(xold1)
        xold1 .= copy(x)
        x .= copy(xmma)

        # Update function values and gradients
        # Note that OPT[:dv] gets updated inside these functions
        f0val , df0dx = obj(FE,OPT,GEOM,x)
        fval = nonlcon(FE,OPT,GEOM,x)[1]
        dfdx = nonlcongrad(FE,OPT,GEOM,x)

        dfdx = transpose(dfdx)      # m x n matrix to pass as dfdx into MMA

        # Compute change in design variables
        # Check only after first iteration
        if iter > 1
            dv_step_change = norm(x .- xold1)
            if dv_step_change < OPT[:options][:step_tol]
                @printf("Design step convergence tolerance satisfied.\n")
                println(fid, "Design step convergence tolerance satisfied.")
            end
        end

        if iter == OPT[:options][:max_iter]
            @printf("Reached maximum number of iterations.\n")
            println(fid, "Reached maximum number of iterations.")
        end

        # Compute norm of KKT residual vector
        residu, kktnorm, residumax = kktcheck(ncons, ndv, xmma, ymma, zmma,
            lam, xsi, eta, mu, zet, s, lb, ub, df0dx, fval, dfdx, a0, a, c, d)

        println(fid, "It. ", iter, ", Obj = ", @sprintf("%.5e", f0val), ", ConsViol = ", @sprintf("%.5e", fval[1]), 
                    ", KKT-norm = ", @sprintf("%.5e", kktnorm), ", DV norm change = ", @sprintf("%.5e", dv_step_change))
        
        # Produce output to screen
        @printf("It. %d, Obj = %.5e, ConsViol = %.5e, KKT-norm = %.5e, DV norm change = %.5e\n",
                iter, f0val, fval[1], kktnorm, dv_step_change)

        push!(history[:fval], f0val)
        push!(history[:fconsval], fval)
        push!(history[:x], x)

        # optional VTK output
        if OPT[:options][:write_to_vtk] == "all"
            writevtk(OPT[:options][:vtk_output_path], "dens", iter)
        end

        plotfun(iter)
    end

    # optional VTK output
    if OPT[:options][:write_to_vtk] == "last"
        writevtk(OPT[:options][:vtk_output_path], "dens", iter)
    end

    println(fid, "The final compliance is : ", round(OPT[:compliance], digits=3))
    println(fid, "The final volume fraction is: ", round(OPT[:volume_fraction], digits=3))
    println(fid, "End time: ", Dates.format(now(), "yyyy-mm-dd E HH:MM:SS"))
    close(fid)

    return history
end