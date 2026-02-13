# =========================================================================
# 
# GPTO Julia with Gridap
# 
# A Julia/Gridap adaptation of the MATLAB code for topology optimization 
# with bars using the geometry projection method, modified for solving
# heat conduction (or heat diffusion) problem.
#
# Version 0.1.0 -- August 2025
#
# This code is a migration to Julia/Gridap by
#
# Prabin Pradhananga
# School of Mechanical, Aerospace, and Manufacturing Engineering
# University of Connecticut
#
# of PyGPTO (Python/Numpy adaptation of the MATLAB code) for topology optimization 
# with bars using the geometry projection method, translated to Python by
#
# Andres Ortegon
# Department of Mathematics
# National Unversity of Colombia
#
# from the GPTO Matlab code written by
#
# Hollis Smith and Julian Norato
# Department of Mechanical Engineering
# University of Connecticut
#
# Special thanks to Sy Nguyen-Van for his Julia code that helped a lot during
# this translation.
#
#
# Disclaimer
# ==========
# This software is provided by the contributors "as-is" with no explicit or
# implied warranty of any kind. In no event shall the University of
# Connecticut or the contributors be held liable for damages incurred by
# the use of this software.
#
# License
# =======
# This software is released under the Creative Commons CC BY-NC 4.0
# license. As such, you are allowed to copy and redistribute the material 
# in any medium or format, and to remix, transform, and build upon the 
# material, as long as you: 
# a) give appropriate credit, provide a link to the license, and indicate 
# if changes were made. You may do so in any reasonable manner, but not in 
# any way that suggests the licensor endorses you or your use.
# b) do not use it for commercial purposes.
#
# To fulfill part a) above, we kindly ask that you please cite the paper
# that introduces this code:
#
# Smith, H. and Norato, J.A. "A MATLAB code for topology optimization
# using the geometry projection method."
# Structural and Multidisciplinary Optimization, 2020,
# https://doi.org/10.1007/s00158-020-02552-0
#
# =========================================================================

## source folders containing scripts not in this folder

using Dates
using Printf, DelimitedFiles
using LinearAlgebra, SparseArrays, StaticArrays
using Makie, GLMakie
using Gridap, Gridap.Geometry, Gridap.Fields, GridapGmsh, Gridap.CellData, Gridap.FESpaces
using Zygote

include("get_inputs.jl")    # get inputs from models
include("FE_routines.jl") 
include("mesh_utilities.jl")
include("geometry_projection.jl")
include("optimization.jl")
include("functions.jl")
include("utilities.jl")
include("MMA.jl")
include("plotting.jl")

## Initialization
init_FE(FE,OPT,GEOM)
init_geometry(FE,OPT,GEOM)
init_optimization(FE,OPT,GEOM)

## Analysis
perform_analysis(FE,OPT,GEOM)

## Finite difference check of sensitivities
# (If requested)
if OPT[:make_fd_check] == true
    run_finite_difference_check(FE,OPT,GEOM)
    exit()   # End code here
end

## Optimization
OPT[:history] = runmma(FE,OPT,GEOM,copy(OPT[:dv]),obj,nonlcon)

## Plot History
if OPT[:options][:plot] == true
    plot_history()
end


