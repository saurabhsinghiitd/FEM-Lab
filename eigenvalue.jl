using LinearAlgebra
using Plots

# Function to plot eigenvectors
function plot_eigenvectors(eigvecs::AbstractMatrix, title_text="Eigenvectors")
    n = size(eigvecs, 1)
    # Defining Plot details
    p = plot(title=title_text, legend=:topright,aspect_ratio=:equal)

    # Add a circle of unit radius on the plot
    plot!(cos, sin, 0, 2π, label="",framestyle=:origin) 
    
    # Plot both eigenvectors on same plot
    for i in 1:n
        plot!(p, [0, eigvecs[1, i]], [0, eigvecs[2, i]] , arrow=true, label="Eigenvector $i")
    end
    
    xlims!(-1, 1)
    ylims!(-1, 1)
    savefig("eigenvectors.png")

end

function main()
    K = [1 2 ; 1 -1] #2D matrix
    eigvals, eigvecs = eigen(K)     #function defined in LinearAlgebra package
    println("Matrix: ")
    display(K)
    println("\nEigenvalues: ")
    display(eigvals)
    println("\nEigenvectors: ")
    display(eigvecs)
    # Plot eigenvectors
    plot_eigenvectors(eigvecs, "Eigenvectors of K")
end

main()