# traffic_dynstab
Traffic Dynamic Stability Research

TODO:
 - Go through and change all matrix names
 - Call generate_grid once and pass the grid as a parameter rather than n
 - Write in print statement to keep track of what the long runs are doing
 - Clean up code if necessary
 - Split the functions into modular files
 - Find where the codes speed is bottlenecking

 - Plot how the maximum eigenvalues increase with n to check the convexity
 - Create interpolation and plotting functions and look at eigenfunctions at some specific n (with satisfactorily low error) to see how they evolve as the eigenvalue increases
 - Write finite difference methods for the general PDE to see what I can find there



IDEAS:
 - Higher order finite difference methods?
 - Binary search to find the nearest grid points?
 - Polynomial interpolation over linear?
 - Recursively define grid points by z^k(x_i) to minimize interpolation screwyness? Nonuniform grid might just make the finite difference stuff worse.
 - choose grid points as Chebyshev points to minimize Runge phenomina when interpolating Y(x)? Maybe I only need to do this when I'm analyzing the eigenvectors. will switching back and forth screw with my results? Will also make finite differences a little worse.
 - Parallelize the operations to reduce runtime? 
 - Penn State offers cloud computing resources, but I think I still have to pay. Not sure what the point of that is then, but maybe the department would pay for it if I ended up really needing more power.

QUESTIONS:
 - Does the finite difference stuff stop working if Y(x) isn't infinitely differentiable?
 - What if the eigenfunction for a given eigenvalue isn't unique? How do we know which eigenfunction we converge to? What would this even mean? Is it even possible? I'm guessing that it's not a problem because our numerics aren't precise enough to overlap eigenfunctions, but analytically I don't see any reason for it to not be possible.
 - We're connsidering the eigenproblem right now, but once I try to solve the general PDE, what happens if there exists no solution for my choice of a(x), b(x), and z(x)? How do the numerical methods interpret that?
 - What if the error terms for the finite difference aren't convergent? What if the nth derivative of Y(x) grows faster than n! and the radius of convergence is smaller than the inverse of the number of dimensions we can consider? This one is more curiousity than practical.

 - If we have some bound for the error on our maximum eigenvalue at some given precision, then can we assert that the maximum eigenvalue is positive if we know that it is greater than the error of the algorithm? Does this question even make sense?

NOTES:
 - I actually do not know that any of this is correct at all. I hope it is. I have no reason to believe it isn't. But it could very well just be spitting out nonsense. I have not checked it thoroughly.