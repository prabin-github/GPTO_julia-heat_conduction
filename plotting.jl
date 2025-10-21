function cross_2_matrix(x, y)
    z = similar(x)
    for i in axes(x, 2)
        vx = SVector{3}(@view(x[:, i]))
        vy = SVector{3}(@view(y[:, i]))
        vz = cross(vx, vy)
        z[:, i] .= vz
    end
    return z
end


function plot_density(iter)
    fig = Figure()
    ax = Axis(fig[1, 1],
            title = "density, $(OPT[:functions][:objective]) = $(@sprintf("%.3f", OPT[:functions][:f][1][:value]))", 
            xlabel="X", ylabel="Y", aspect=DataAspect())
    image!(ax, reshape(OPT[:penalized_elem_dens], Tuple(FE[:mesh_input][:elements_per_side])), colormap=:jet)
    Makie.resize_to_layout!(fig)
    display(fig)

    # Save design plots
    Makie.save(joinpath(OPT[:options][:plot_output_path], OPT[:problem] * "_" * OPT[:functions][:objective] * "_" * 
                OPT[:functions][:constraints] * "_density_$(iter).png"), fig)
    sleep(0.1)
end


function plot_density_levelsets()
    # global FE, OPT, GEOM

    if FE[:mesh_input][:type] != "generate" && FE[:mesh_input][:type] != "read-home-made"
        println("not yet implemented for non meshgrid conforming meshes")
        return
    end

    ## Change here whether you want to plot the penalized (i.e., effective) or 
    ## the unpenalized (i.e., projected) densities.  By default, we plot the 
    ## effective densities.
    #
    # For penalized, use OPT.penalized_elem_dens;
    # For unpenalized, use OPT.elem_dens;

    plot_dens = OPT[:elem_dens]
    # plot_dens = OPT[:penalized_elem_dens]

    if FE[:dim] == 2
        # Level sets
        n = 64
        levels = range(0, stop=1, length=n)

        if !haskey(OPT[:options], :centroid_mesh)
            OPT[:options][:centroid_mesh] = Dict{Symbol, Any}()

            mn = FE[:mesh_input][:elements_per_side]
            nm = mn[[1,2]]  # for meshgrid, Julia is 1-based

            OPT[:options][:centroid_mesh][:shape] = nm
            OPT[:options][:centroid_mesh][:X] = reshape(FE[:centroids][1, :], Tuple(nm))
            OPT[:options][:centroid_mesh][:Y] = reshape(FE[:centroids][2, :], Tuple(nm))
        end

        X = OPT[:options][:centroid_mesh][:X]
        Y = OPT[:options][:centroid_mesh][:Y]
        V = reshape(plot_dens, Tuple(OPT[:options][:centroid_mesh][:shape]))

        fig = Makie.Figure()
        ax = Axis(fig[1, 1],
                title = "density, $(OPT[:functions][:objective]) = $(@sprintf("%.3f", OPT[:functions][:f][1][:value]))", 
                xlabel="X", ylabel="Y")
        Makie.contourf!(ax, X, Y, V; levels=collect(levels), colormap=Reverse(:grays))
        Makie.resize_to_layout!(fig)
        xlims!(ax, FE[:coord_min][1], FE[:coord_max][1])
        ylims!(ax, FE[:coord_min][2], FE[:coord_max][2])
        display(fig)
        sleep(0.1)

    elseif FE[:dim] == 3
        levels = [0.25, 0.5, 0.75]
        # 3D plotting not implemented here
    end
end


function plot_design(args...)
    # Plot_design(fig,point_mat,bar_mat) plots the bars into the figure fig
    # fig is the number (or handle) of the figure to use
    # global GEOM, FE

    nargs = length(args)

    if nargs == 0
        iter = 1
        point_mat = GEOM[:current_design][:point_matrix]
        bar_mat = GEOM[:current_design][:bar_matrix]
    elseif nargs == 1
        iter = args[1]
        point_mat = GEOM[:current_design][:point_matrix]
        bar_mat = GEOM[:current_design][:bar_matrix]
    elseif nargs == 3
        iter = args[1]
        point_mat = args[2]
        bar_mat = args[3]
    else
        println("plot_design received an invalid number of arguments.")
        return
    end

    ## user specified parameters

    # set the color of the bars
    bar_color = [1 0 0]    # red 
    # set size variable threshold to plot bars
    size_tol = 0.05
    # set the resolution of the bar-mesh (>=8 and even)
    N = 16

    ## bar points,vectors and length
    bar_tol = 1e-12 # threshold below which bar is just a circle
    n_bar = size(bar_mat, 1)

    x_1b = zeros(3, n_bar)
    x_2b = zeros(3, n_bar) # these are always in 3D 

    pt1_IDs = bar_mat[:, 2]
    pt2_IDs = bar_mat[:, 3]

    x_1b[1:FE[:dim], :] = transpose(point_mat[Int.(vec(GEOM[:point_mat_row][pt1_IDs])), 2:end])
    x_2b[1:FE[:dim], :] = transpose(point_mat[Int.(vec(GEOM[:point_mat_row][pt2_IDs])), 2:end])

    n_b = x_2b - x_1b
    l_b = sqrt.(sum(n_b .* n_b, dims=1))  # length of the bars  # (1, n_bar)

    ## principle bar direction
    e_hat_1b = n_b ./ l_b
    short = l_b .< bar_tol
    if any(short)
        e_hat_1b[:, vec(short)] = repeat([1; 0; 0], 1, sum(short))
    end
        
    # determine coordinate direction most orthogonal to bar
    case_1 = abs.(n_b[1, :]) .< abs.(n_b[2, :]) .&& abs.(n_b[1, :]) .< abs.(n_b[3, :])
    case_2 = abs.(n_b[2, :]) .< abs.(n_b[1, :]) .&& abs.(n_b[2, :]) .< abs.(n_b[3, :])
    case_3 = .~(case_1 .| case_2)

    ## secondary bar direction
    e_alpha = zeros(size(n_b))
    e_alpha[1, case_1] .= 1
    e_alpha[2, case_2] .= 1
    e_alpha[3, case_3] .= 1

    e_2b = l_b .* cross_2_matrix(e_alpha, e_hat_1b)
    norm_e_2b = sqrt.(sum(e_2b .^ 2, dims=1))
    e_hat_2b = e_2b ./ norm_e_2b

    ## tertiary bar direction
    e_3b = cross_2_matrix(e_hat_1b, e_hat_2b)
    norm_e_3b = sqrt.(sum(e_3b .^ 2, dims=1))
    e_hat_3b = e_3b ./ norm_e_3b

    ## Jacobian transformation (rotation) matrix R
    R_b = zeros(3, 3, n_bar)
    R_b[:, 1, :] = e_hat_1b
    R_b[:, 2, :] = e_hat_2b
    R_b[:, 3, :] = e_hat_3b

    ## create the reference-sphere mesh
    # if FE.dim == 3
    #     [x,y,z] = sphere(N)
    #     sx1 = z[1:N/2,:]
    #     sy1 = x[1:N/2,:]
    #     sz1 = y[1:N/2,:]
    #     sx2 = z[N/2+1:end,:]
    #     sy2 = x[N/2+1:end,:]
    #     sz2 = y[N/2+1:end,:]
    #     X1 = [sx1[:], sy1[:], sz1[:]]'
    #     X2 = [sx2[:], sy2[:], sz2[:]]'
    # else
    N = N^2
    t = range(-π / 2, stop=-π / 2 + 2π, length=N + 1) |> collect
    x = -cos.(t)
    y = sin.(t)
    z = zeros(size(t))

    cxo = x[Int.(1:N/2)]
    cyo = y[Int.(1:N/2)]
    czo = z[Int.(1:N/2)]

    cxf = x[Int.(N/2+1:end)]
    cyf = y[Int.(N/2+1:end)]
    czf = z[Int.(N/2+1:end)]

    X1 = vcat(cxo[:]', cyo[:]', czo[:]')
    X2 = vcat(cxf[:]', cyf[:]', czf[:]')
    # end

    ## create the surface for each bar and plot it
    r_b = bar_mat[:, end]
    alpha = bar_mat[:, end-1]

    Color = bar_color

    fig = Figure()
    ax = Axis(fig[1, 1], title="design, iteration = " * string(iter), xlabel="X", ylabel="Y", aspect=DataAspect())

    for b = 1:n_bar
        Alpha = alpha[b] .^ 2

        bar_X1 = r_b[b] * R_b[:, :, b] * X1 .+ x_1b[:, b]
        bar_X2 = r_b[b] * R_b[:, :, b] * X2 .+ x_2b[:, b]
        # if FE.dim == 3
        #     bar_x1 = reshape(bar_X1[1,:], [N/2, N+1])
        #     bar_y1 = reshape(bar_X1[2,:], [N/2, N+1])
        #     bar_z1 = reshape(bar_X1[3,:], [N/2, N+1])

        #     bar_x2 = reshape(bar_X2[1,:], [N/2+1, N+1])
        #     bar_y2 = reshape(bar_X2[2,:], [N/2+1, N+1])
        #     bar_z2 = reshape(bar_X2[3,:], [N/2+1, N+1]) 
        # else
        bar_x1 = bar_X1[1, :]
        bar_y1 = bar_X1[2, :]
        bar_z1 = bar_X1[3, :]

        bar_x2 = bar_X2[1, :]
        bar_y2 = bar_X2[2, :]
        bar_z2 = bar_X2[3, :]
        # end

        bar_x = [bar_x1; bar_x2]
        bar_y = [bar_y1; bar_y2]
        bar_z = [bar_z1; bar_z2]

        if Alpha > size_tol
            # if FE.dim == 3
            #     s = surfl(bar_x,bar_y,bar_z) # shaded surface with lighting
            #     s.LineStyle = 'none'
            #     s.FaceAlpha = Alpha
            #     shading interp
            # else
            # s = patch(bar_x,bar_y,colormap =:jet)
            # s.FaceAlpha = Alpha
            # Plot the polygon using bar_x and bar_y
            poly!(ax, bar_x, bar_y, color=(:red, Alpha), strokecolor=:black, strokewidth=1.5, overdraw=true)
            # end
        end
    end
    Makie.resize_to_layout!(fig)
    xlims!(ax, FE[:coord_min][1], FE[:coord_max][1])
    ylims!(ax, FE[:coord_min][2], FE[:coord_max][2])
    display(fig)

    # Save design plots
    Makie.save(joinpath(OPT[:options][:plot_output_path], OPT[:problem] * "_" * OPT[:functions][:objective] * "_" * 
                OPT[:functions][:constraints] * "_design_$(iter).png"), fig)
    sleep(0.1)
end


function plot_history()
    # Top subplot: Objective history (semilogy)
    fig = Figure();
    ax1 = Axis(fig[1, 1], title="Objective and Constraint History", yscale=log10)

    # Pull out the vector of objective values – assuming stored as a 1×N or N×1 array or Vector
    fvals = vec(OPT[:history][:fval])
    lines!(ax1, 0:length(fvals)-1, fvals, label=OPT[:functions][:f][1][:name])
    axislegend(ax1, position=:rt)

    # Bottom subplot: Constraint history (if present)
    if haskey(OPT[:history], :fconsval)
        # reshape into (n_iter × n_constr)
        g = reshape(vec(OPT[:history][:fconsval]), :, OPT[:functions][:n_func]-1) .+ 
                OPT[:functions][:constraint_limit]

        # build labels & scales
        labels = String[]
        scales = ones(OPT[:functions][:n_func]-1)
        for i in 2:OPT[:functions][:n_func]
            push!(labels, OPT[:functions][:f][i][:name])
            if OPT[:functions][:f][i][:name] == "angle constraint"
                scales[i-1] = OPT[:options][:angle_constraint][:scale]
            end
        end

        ax2 = Axis(fig[2, 1]; xlabel="Iteration")
        # Plot each constraint curve
        for i in axes(g,2)
            lines!(ax2, 0:size(g,1)-1, g[:,i]./scales[i], label=labels[i], color=:cyan)
        end

        # Dashed horizontal line at the limit
        # hlines!(ax2, OPT[:functions][:constraint_limit], 
        #         linestyle=:dash, color=:black, label="limit")
        lines!(ax2,
            [0, size(g,1)-1],           # x‐coordinates
            [OPT[:functions][:constraint_limit][1], OPT[:functions][:constraint_limit][1]],     # y‐coordinates
            linestyle=:dash,
            color=:black,
            label="limit")
        axislegend(ax2, position=:rb)
    end

    # Save iteration history plots
    Makie.save(joinpath(OPT[:options][:vtk_output_path], OPT[:problem] * "_" * OPT[:functions][:objective] * "_" * 
                OPT[:functions][:constraints] * "_History.png"), fig)
    
    display(GLMakie.Screen(), fig)
    sleep(2)
end


function writevtk(folder, name_prefix, iteration)
    # This function writes a vtk file with the mesh and the densities that can
    # be plotted with, e.g., ParaView
    #
    # This function writes an unstructured grid (vtk format) to folder (note
    # that the folder is relative to the rooth folder where the main script is
    # located).
    #
    # NOTE: if a vtk file with the same name exists in the folder, it will be
    # overwritten.
    # global FE, OPT

    # Make sure the output folder exists, and if not, create it
    if !isdir(folder)
        mkpath(folder)
    end

    # build filename
    suffix   = @sprintf("%03d", iteration)
    filename = joinpath(folder, string(name_prefix, suffix, ".vtk"))

    # Open the file for writing (overwrite if exists)
    open(filename, "w") do fid
        # Write header
        println(fid, "# vtk DataFile Version 1.0")
        println(fid, "Bar_TO_3D")
        println(fid, "ASCII")
        println(fid, "DATASET UNSTRUCTURED_GRID")

        # Write nodal coordinates
        coords = zeros(3, FE[:n_node])
        coords[1:FE[:dim], :] = FE[:coords][1:FE[:dim], :]

        println(fid, "POINTS " * string(FE[:n_node]) * " float ")
        for inode in 1:FE[:n_node]
                println(fid, coords[1, inode], " ", coords[2, inode], " ", coords[3, inode], " ")
        end

        # Write elements
        nnodes = 2^FE[:dim]  # 4 for quads, 8 for hexas

        println(fid, " CELLS " * string(FE[:n_elem]) * " " * string(FE[:n_elem] * (nnodes + 1)) * "  ")
        
        for iel in 1:FE[:n_elem]
            # IMPORTANT! Vtk numbers nodes from 0, so we subtract 1
            if FE[:dim] == 2
                nel = 4
            elseif FE[:dim] == 3
                nel = 8
            end
            writedlm(fid, [nel vec(FE[:elem_node][:, iel])' .- 1], " ")
        end

        # Write element types
        println(fid, "CELL_TYPES " * string(FE[:n_elem]) * " ")

        if FE[:dim] == 2
            elem_type = 9  # Corresponding to VTK_QUAD
        elseif FE[:dim] == 3
            elem_type = 12 # Corresponding to VTK_HEXAHEDRON
        end

        for iel in 1:FE[:n_elem]
            println(fid, elem_type)
        end
        
        # Write elemental densities
        println(fid, "CELL_DATA " * string(FE[:n_elem]) * " ")
        println(fid, "SCALARS density float 1 ")
        println(fid, "LOOKUP_TABLE default ")
        for iel in 1:FE[:n_elem]
            density = OPT[:elem_dens][iel]
            println(fid, density, " ")
        end
    end

end
