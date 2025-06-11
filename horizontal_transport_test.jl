# Example verification of horizontal transport coefficient placement
#
# This script constructs a 2x2 grid (two altitude layers, two horizontal columns)
# with a single species transported horizontally.  Chemistry is disabled so the
# only processes are horizontal advection.  The transport matrix is built using
# the same coefficient placement as bluejay's `chemJmat` and `ratefn` functions.
#
# The resulting rate matrix can be compared against the analytic solution
# produced in `coupled_cells_analytic_solution.nb`.

# using LinearAlgebra

# # Grid setup ---------------------------------------------------------------
# num_alt = 2      # two altitude bins
# num_horiz = 2    # two horizontal columns
# Δx = 1.0         # cm, horizontal cell width

# # Constant wind blowing from column 1 toward column 2 at all altitudes
# v = 10.0         # cm/s

# # Transport coefficients follow update_horiz_transport_coefficients()
# adv_rate = v / Δx

# # Transport matrix for one altitude layer
# T = [-adv_rate  adv_rate;
#       adv_rate -adv_rate]

# # Build full 4×4 matrix for two independent altitude layers
# A = kron(Matrix(I, num_alt, num_alt), T)

# # Example initial state: higher density in column 1
# n₀ = [1.0, 0.0, 1.0, 0.0]

# # Time derivative using the transport matrix
# println("d n/dt at t=0:")
# println(A * n₀)

# # Jacobian equals the transport matrix for this linear problem
# println("Jacobian matrix:")
# println(A)

# # Analytic solution after a short time step Δt
# Δt = 0.1
# n_t = exp(Δt * A) * n₀
# println("n(t=Δt):")
# println(n_t)

###########################################################################################################################################

# Example verification of horizontal transport coefficient placement
#
# This script constructs a 2x2 grid (two altitude layers, two horizontal columns)
# with a single species transported horizontally.  Chemistry is disabled so the
# only processes are horizontal advection.  The transport matrix is built using
# the same coefficient placement as bluejay's `chemJmat` and `ratefn` functions.
#
# The resulting rate matrix can be compared against the analytic solution
# produced in `coupled_cells_analytic_solution.nb`.

# using LinearAlgebra

# # Grid setup ---------------------------------------------------------------
# num_alt = 2      # two altitude bins
# num_horiz = 2    # two horizontal columns
# Δx = 1.0         # cm, horizontal cell width

# # Constant wind blowing from column 1 toward column 2 at all altitudes
# v = 10.0         # cm/s

# # Transport coefficients follow update_horiz_transport_coefficients()
# adv_rate = v / Δx

# # Build transport matrix using the same indexing scheme as chemJmat
# function build_matrix(num_alt, num_horiz, adv)
#     n = num_alt * num_horiz
#     A = zeros(n, n)
#     entries = Tuple{Int,Int,Float64}[]
#     for ih in 1:num_horiz
#         for ia in 1:num_alt
#             idx = (ih - 1) * num_alt + ia
#             if ih != num_horiz
#                 j = ih * num_alt + ia
#                 A[idx, idx] -= adv
#                 A[idx, j] += adv
#                 push!(entries, (idx, idx, -adv))
#                 push!(entries, (idx, j, adv))
#             end
#         end
#     end
#     return A, entries
# end

# A, entries = build_matrix(num_alt, num_horiz, adv_rate)

# # Analytic matrix expected from the Mathematica solution
# A_expected = [-adv_rate 0 adv_rate 0;
#               0 -adv_rate 0 adv_rate;
#               adv_rate 0 -adv_rate 0;
#               0 adv_rate 0 -adv_rate]

# # Example initial state: higher density in column 1
# n₀ = [1.0, 0.0, 1.0, 0.0]

# # Time derivative using the transport matrix
# println("d n/dt at t=0:")
# println(A * n₀)

# # Jacobian equals the transport matrix for this linear problem
# println("Jacobian matrix:")
# println(A)
# println("Matches analytic matrix: ", A ≈ A_expected)

# println("Sparse triplets (i,j,value):")
# for (i,j,v) in entries
#     println("(", i, ", ", j, ", ", v, ")")
# end

# # Analytic solution after a short time step Δt
# Δt = 0.1
# n_t = exp(Δt * A) * n₀
# println("n(t=Δt):")
# println(n_t)

###########################################################################################################################################

# Example verification of horizontal transport coefficient placement
#
# This script constructs a 2x2 grid (two altitude layers, two horizontal columns)
# with a single species transported horizontally.  Chemistry is disabled so the
# only processes are horizontal advection.  The transport matrix is built using
# the same coefficient placement as bluejay's `chemJmat` and `ratefn` functions.
#
# The resulting rate matrix can be compared against the analytic solution
# produced in `coupled_cells_analytic_solution.nb`.

using LinearAlgebra

# Grid setup ---------------------------------------------------------------
num_alt = 2      # two altitude bins
num_horiz = 2    # two horizontal columns
Δx = 1.0         # cm, horizontal cell width

# Constant wind blowing from column 1 toward column 2 at all altitudes
v = 10.0         # cm/s

# Transport coefficients follow update_horiz_transport_coefficients()
adv_rate = v / Δx

# Build transport matrix using the same indexing scheme as chemJmat
# In bluejay, altitude varies more slowly than the horizontal index when
# densities are flattened. The index for altitude `ia` and column `ih` is
# `(ia - 1) * num_horiz + ih`.
function build_matrix(num_alt, num_horiz, adv)
    n = num_alt * num_horiz
    A = zeros(n, n)
    entries = Tuple{Int,Int,Float64}[]

    for ia in 1:num_alt           # loop altitude first
        for ih in 1:num_horiz     # columns vary fastest
            idx = (ia - 1) * num_horiz + ih
            if ih < num_horiz
                idx_right = idx + 1
                A[idx, idx] -= adv
                A[idx, idx_right] += adv
                push!(entries, (idx, idx, -adv))
                push!(entries, (idx, idx_right, adv))

                # symmetric entry for the cell to the right
                A[idx_right, idx_right] -= adv
                A[idx_right, idx] += adv
                push!(entries, (idx_right, idx_right, -adv))
                push!(entries, (idx_right, idx, adv))
            end
        end
    end

    return A, entries
end

A, entries = build_matrix(num_alt, num_horiz, adv_rate)

# Analytic matrix expected from the Mathematica solution
A_expected = [-adv_rate 0 adv_rate 0;
              0 -adv_rate 0 adv_rate;
              adv_rate 0 -adv_rate 0;
              0 adv_rate 0 -adv_rate]

# Example initial state: higher density in column 1
n₀ = [1.0, 0.0, 1.0, 0.0]

# Time derivative using the transport matrix
println("d n/dt at t=0:")
println(A * n₀)

# Jacobian equals the transport matrix for this linear problem
println("Jacobian matrix:")
println(A)
println("Matches analytic matrix: ", A ≈ A_expected)

println("Sparse triplets (i,j,value):")
for (i,j,v) in entries
    println("(", i, ", ", j, ", ", v, ")")
end

# Analytic solution after a short time step Δt
Δt = 0.1
n_t = exp(Δt * A) * n₀
println("n(t=Δt):")
println(n_t)

