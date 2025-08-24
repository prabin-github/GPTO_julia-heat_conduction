function fd_check_cost(FE,OPT,GEOM)
    #
    # This function performs a finite difference check of the sensitivities of
    # the COST function with respect to the bar design variables.
    # global OPT 
    #
    # ===============================
    # FINITE DIFFERENCE SENSITIVITIES
    # ===============================
    n_dv = OPT[:n_dv]
    grad_theta_i = zeros(n_dv, 1)

    fd_step = OPT[:fd_step_size]

    max_error = 0.0
    max_rel_error = 0.0
    max_error_dv = 0
    max_rel_error_dv = 0

    dv_0 = copy(OPT[:dv])
    dv_i = copy(OPT[:dv])

    theta_0, grad_theta_0 = obj(FE,OPT,GEOM,dv_0)

    # Finite differences
    println("Computing finite difference sensitivities of cost function...")

    # Do this for all design variables or only a few
    up_to_dv = n_dv

    for i in 1:up_to_dv
        # Preturb dv
        dv_i[i] = dv_0[i] + fd_step
        theta_i, _ = obj(FE,OPT,GEOM,dv_i)
        grad_theta_i[i] = (theta_i - theta_0) / fd_step
        error = grad_theta_0[i] - grad_theta_i[i]

        if abs(error) > abs(max_error)
            max_error = error
            max_error_dv = i
        end
        rel_error = error / theta_0
        if abs(rel_error) > abs(max_rel_error)
            max_rel_error = rel_error
            max_rel_error_dv = i
        end
        dv_i = copy(dv_0)
    end

    # theta_0, grad_theta_0 = obj(FE,OPT,GEOM,dv_0)  # to reset the design
    OPT[:dv] = dv_0

    @printf("Max. ABSOLUTE error is: %.5e\n", max_error)
    @printf("It occurs at variable: %d\n", max_error_dv)

    @printf("Max. RELATIVE error is: %.5e\n", max_rel_error)
    @printf("It occurs at variable: %d\n", max_rel_error_dv)

    n_dv_vec = LinRange(1, OPT[:n_dv], OPT[:n_dv])  # Vector of index of dv
    name_obj = OPT[:problem] * ": FD of cost function: " * OPT[:functions][:objective]  # Name of plots
    fig_FD_cost = Figure()
    ax_FD_cost = Axis(fig_FD_cost[1, 1], xlabel="Design variable: v", ylabel="dz/dv", title=name_obj)
    Makie.lines!(ax_FD_cost, n_dv_vec, vec(grad_theta_0), color=:red, label="Analytical")
    Makie.scatter!(ax_FD_cost, n_dv_vec, vec(grad_theta_i), color=:blue, label="FD")
    axislegend(position=:rb)
    display(fig_FD_cost)
    Makie.save(OPT[:problem] * "_FD_cost_" * OPT[:functions][:objective] * ".png", fig_FD_cost)  # Saving plots
end


function fd_check_constraint(FE,OPT,GEOM)
    # !!!!!!!!!!!!!!!!!! need to work for number of constraints n_con
    #
    # This function performs a finite difference check of the sensitivities of
    # the CONSTRAINT function with respect to the bar design variables.
    # It is currently setup for one constraint, but it can be easily modified
    # for other/more constraints.
    # global OPT

    # ===============================
    # FINITE DIFFERENCE SENSITIVITIES
    # ===============================
    n_dv = OPT[:n_dv]
    grad_theta_i = zeros(n_dv, 1)

    fd_step = OPT[:fd_step_size]

    max_error = 0.0
    max_rel_error = 0.0
    max_error_dv = 0
    max_rel_error_dv = 0

    dv_0 = copy(OPT[:dv])
    dv_i = copy(OPT[:dv])

    theta_0 = nonlcon(FE,OPT,GEOM,dv_0)
    grad_theta_0 = nonlcongrad(FE,OPT,GEOM,dv_0)

    # Finite differences
    println("Computing finite difference sensitivities of constraint function...")
    
    # Do this for all design variables or only a few
    up_to_dv = n_dv

    for i in 1:up_to_dv
        #perturb dv
        dv_i[i] = dv_0[i] + fd_step
        theta_i = nonlcon(FE,OPT,GEOM,dv_i)
       
        grad_theta_i[i] = (theta_i[1] .- theta_0[1]) ./ fd_step

        error = grad_theta_0[i] - grad_theta_i[i]
       
        if abs(error) > abs(max_error)
            max_error = error
            max_error_dv = i
        end

        rel_error = error / theta_0[1]

        if abs(rel_error) > abs(max_rel_error)
            max_rel_error = rel_error
            max_rel_error_dv = i
        end

        dv_i = copy(dv_0)
    end

    OPT[:dv] = dv_0

    @printf("Max. ABSOLUTE error is: %.5e\n", max_error)
    @printf("It occurs at variable:   %d\n", max_error_dv)

    @printf("Max. RELATIVE error is: %.5e\n", max_rel_error)
    @printf("It occurs at variable: %d\n", max_rel_error_dv)

    n_dv_vec = LinRange(1, OPT[:n_dv], OPT[:n_dv])  # Vector of index of dv
    name_const = OPT[:problem] * ": FD of constraint function: " * OPT[:functions][:constraints]  # Name of plots
    fig_FD_const = Figure()
    ax_FD_const = Axis(fig_FD_const[1, 1], xlabel="Design variable: v", ylabel="dz/dv", title=name_const)
    Makie.lines!(ax_FD_const, n_dv_vec, vec(grad_theta_0), color=:red, label="Analytical")
    Makie.scatter!(ax_FD_const, n_dv_vec, vec(grad_theta_i), color=:blue, label="FD")
    axislegend(position=:rb)
    display(fig_FD_const)
    Makie.save(OPT[:problem] * "_FD_const_" * OPT[:functions][:constraints] * ".png", fig_FD_const)  # Saving plots
end


function run_finite_difference_check(FE,OPT,GEOM)
    # This function performs a finite difference check of the analytical
    # sensitivities of the cost and/or constraint functions by invoking the
    # corresponding routines.
    # global OPT
    if OPT[:check_cost_sens] == true
        fd_check_cost(FE,OPT,GEOM)
    end
    if OPT[:check_cons_sens] == true
        fd_check_constraint(FE,OPT,GEOM)
    end
end