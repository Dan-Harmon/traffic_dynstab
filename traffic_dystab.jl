using SparseArrays, Arpack, LinearMaps

function generate_grid(n::Int; xl=0, xh=1)::Tuple{Array{Float64}, Float64}
#= 
    Generates the discrete grid for us to evaluate the functions Y(x) at the points of. 

    Inputs:
        * n - represents the number of grid points -1
        * xl - the lower bound for the interval upon which we want to construct our grid
        * xh - the upper bound for the interval upon which we want to construct our grid
    Outputs:
        * collect(x) - array of grid values distributed uniformly across the interval [xl,xh]
        * dx - the distance between each grid point
    Ideas:
        * might want to recursively define the grid points by z^k(x_i) to minimize the funky approximations on the off grid points, might make things worse though, now that I'm thinking about it 
=#
    dx = 1/n
    x = xl:dx:1xh

    return collect(x), dx
end

function finite_diff_matrix(n::Int)::SparseMatrixCSC{Float64}
#=
    Constructs an n+1xn+1 dimension finite difference matrix to approximate Y'(x) on the grid points

    Inputs:
        * n - represents the dimension of the finite difference matrix -1
    Outputs:
        * The n+1xn+1 finite differnce matrix used to approximate Y'(x)
=#
    x, dx = generate_grid(n)

    A = spzeros(n+1, n+1)
    
    A[1,1] = -1/dx
    A[1,2] = 1/dx

    A[n+1,n] = -1/dx
    A[n+1,n+1] = 1/dx

    for i in 2:n
        A[i,i-1] = -1/(2*dx)
        A[i,i+1] = 1/(2*dx)
    end

    return A
end

function offgrid_finite_diff_matrix(n::Int, z::Function)::SparseMatrixCSC{Float64}
#=
    Constructs an n+1xn+1 dimension finite difference matrix to approximate Y'(z(x)) at the off grid points

    Inputs:
        * n - represents the dimension of the finite difference matrix -1
        * c - the coeffient in (0,1) defining z(x)=cx
    Outputs:
        * The n+1xn+1 finite differnce matrix used to approximate Y'(z(x))
=#
    x, dx = generate_grid(n)

    B = spzeros(n+1,n+1)

    for i in 1:(n+1)
        tol = 1e-10
        maybe_match = findfirst(y -> abs(z(x[i]) - y) < tol, x)

        if maybe_match !== nothing
            j = maybe_match
            if j == 1
                B[i,1] = -1/dx
                B[i,2] = 1/dx
            elseif j == n+1
                B[i,n] = -1/dx
                B[i,n+1] = 1/dx
            else
                B[i,j-1] = -1/(2*dx)
                B[i,j+1] = 1/(2*dx) 
            end
        else
            idx_below = findlast(y -> y <= z(x[i]), x)
            idx_above = idx_below + 1

            B[i, idx_below] = -1/dx
            B[i, idx_above] = 1/dx
        end
    end

    return B
end

function construct_grid_matrix(n::Int, f::Function)::SparseMatrixCSC{Float64}
#=
    Constructs a diagonal matrix represnting the values of some function f(x_i) at each grid point

    Inputs:
        * n - represents the dimension of the derivative matrix -1
        * f - the function which we are evaluating at the grid points
    Outputs:
        * M - the diagonal matrix constaining the values of f(x_i) at each grid point
=#
    x, dx = generate_grid(n)

    M = spdiagm(0 => f.(x))
    
    return M
end

function find_max_eval(Mmap::LinearMap; nev=1)::Complex{Float64}
#=
    Solves the for the largest eigenvalue for the linear transformation A+cB
        
    Inputs:
        * Mmap - an abstract linear map representing a'(x)Y(x)+a(x)Y'(x)-z'(x)Y'(x) = (DPHA+PHA*A-DZ*B)Y(x) evaluated at the grid points
    Outputs:
        * λ_max - The largest eigenvalue of the linear transformation Mmap
=#

    λ_max = eigs(Mmap, nev=nev, which=:LR, maxiter=1000, tol=1e-6)[1][1]
    return λ_max
end

function make_operator(A::SparseMatrixCSC{Float64}, B::SparseMatrixCSC{Float64}, DZ::SparseMatrixCSC{Float64}, PHA::SparseMatrixCSC{Float64}, DPHA::SparseMatrixCSC{Float64}, ETA::SparseMatrixCSC{Float64}, DETA::SparseMatrixCSC{Float64}, YZ::SparseMatrixCSC{Float64})::LinearMap
#= 
    Function to construct the linear operator DPHA+PHA*A-DZ*B, saves memory from explicitly saving the matrix as n grows large
    Inputs:
        * A - the finite difference matrix approximating Y'(x) at the grid points
        * B - the finite difference matrix approximating Y'(z(x)) at the grid points
        * DZ - the diagonal matrix represting z'(x) at the grid points
        * PHA - the diagonal matrix represting a(x) at the grid points
        * DPHA - the diagonal matrix represting a'(x) at the grid points
        * ETA - the diagonal matrix represting b(x) at the grid points
        * DETA - the diagonal matrix represting b'(x) at the grid points
        * YZ - the linear interpolation matrix  approximating Y(z(x)) at the grid points
    Outputs:
        * LinearMap - an abstract mapping to represent the linear operator DPHA+PHA*A-DZ*B without storing an explicit n+1xn+1 matrix
=#

    return LinearMap(v -> (DPHA + DETA*YZ + PHA*A - ETA*DZ*B)*v, size(A,1), size(A,1))
end

function linear_interpolation(n::Int, f::Function)::SparseMatrixCSC{Float64}
#=
    Generates the linear interpolation for Y(z(x)) from the two nearest grid points 
    
    Inputs:
        * n - represents the dimension of the derivative matrix -1
        * f - the function which we are evaluating at the grid points
    Outputs:
        * YZ - the matrix representing the interpolation Y(z(x)) = YZ * Y(x) on the grid points
=#
    x, dx = generate_grid(n)

    YZ = spzeros(n+1,n+1)

    for i in 1:(n+1)
        tol = 1e-10
        maybe_match = findfirst(y -> abs(z(x[i]) - y) < tol, x)

        if maybe_match !== nothing
            j = maybe_match
            YZ[i,j] = 1   
        else
            idx_below = findlast(y -> y <= z(x[i]), x)
            idx_above = idx_below + 1

            # I am reasonably sure this is correct, but honestly I let myself get lost in the sauce a bit
            YZ[i, idx_below] = (x[idx_above]-z(x[i]))/dx
            YZ[i, idx_above] = (z(x[i])-x[idx_below])/dx
        end
    end

    return YZ   

end



z(x) = 0.5(x+ x^2)
dz(x) = 0.5+x
a(x) = exp(-2*x)
da(x) = -2*exp(-2*x)

b(x) = 1
db(x) = 0



for n in [10, 20, 50, 100, 200, 300, 500, 800, 1200]
    A = finite_diff_matrix(n)
    B = offgrid_finite_diff_matrix(n, z)
    PHA = construct_grid_matrix(n, a)
    DPHA = construct_grid_matrix(n, da)
    ETA = construct_grid_matrix(n, b)
    DETA = construct_grid_matrix(n, db)
    DZ = construct_grid_matrix(n, dz)
    YZ = linear_interpolation(n, z)
    Mmap = make_operator(A,B,DZ, PHA, DPHA, ETA, DETA, YZ)
    λ_max = find_max_eval(Mmap)
    print("n=$n, λ_max = $λ_max \n\n")
end
