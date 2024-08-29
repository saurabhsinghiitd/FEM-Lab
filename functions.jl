module functions
using Plots
using LinearAlgebra
function formStiffness2Dtruss(GDof, numberElements, elementNodes, nodeCoordinates, E_vec, A_vec)
    # Initialize the global stiffness matrix
    K_Assembly = zeros(GDof, GDof)
    for iElement in 1:numberElements
        iNodes = elementNodes[iElement, :]
        elementDof = [iNodes[1] * 2 - 1, iNodes[1] * 2, iNodes[2] * 2 - 1, iNodes[2] * 2]
  
        D_x = nodeCoordinates[iNodes[2], 1] - nodeCoordinates[iNodes[1], 1]
        D_y = nodeCoordinates[iNodes[2], 2] - nodeCoordinates[iNodes[1], 2]
        L = sqrt(D_x^2 + D_y^2)
        l = D_x / L
        m = D_y / L
  
        k_global = (E_vec[iElement] * A_vec[iElement] / L) * [
            l^2    l*m    -l^2   -l*m;
            l*m    m^2    -l*m   -m^2;
            -l^2   -l*m    l^2    l*m;
            -l*m   -m^2    l*m    m^2
        ]
  
        K_Assembly[elementDof, elementDof] += k_global
    end
    return K_Assembly
end
  
function solution(prescribedDof, K_assembly, D_vec, F_vec, F_eq_vec=zeros(size(F_vec)))
    GDof = length(D_vec)
  
    freeDof = setdiff(1:GDof, prescribedDof)
  
    # Solve for the displacement vector at free degrees of freedom
    D_vec[freeDof] = K_assembly[freeDof, freeDof] \ (F_vec[freeDof] + F_eq_vec[freeDof])
  
    # Find non-zero degrees of freedom
    nonZeroDof = union(freeDof, prescribedDof[D_vec[prescribedDof] .!= 0])
  
    # Compute the force vector
    F_vec[prescribedDof] = K_assembly[prescribedDof, nonZeroDof] * D_vec[nonZeroDof] - F_eq_vec[prescribedDof]
    return D_vec, F_vec
end
  
function stresses2Dtruss(N_elements, elementNodes, nodesCoordinates, D_col, E_vec)
    stress = fill(NaN, N_elements)
    for iElement in 1:N_elements
        iNodes = elementNodes[iElement, :]
        elementDofs = [iNodes[1] * 2 - 1, iNodes[1] * 2, iNodes[2] * 2 - 1, iNodes[2] * 2]
  
        xa = nodesCoordinates[iNodes[2], 1] - nodesCoordinates[iNodes[1], 1]
        ya = nodesCoordinates[iNodes[2], 2] - nodesCoordinates[iNodes[1], 2]
        L = sqrt(xa^2 + ya^2)
        C = xa / L
        S = ya / L
  
        stress[iElement] = (E_vec[iElement] / L) * dot([-C -S C S],D_col[elementDofs])
    end
    return stress
end
  
  
# Function to plot the mesh
function drawingMesh(plot1, nodeCoords, elementNodes, lc,ls)
    for iElement in 1:size(elementNodes, 1)
        xCoords = nodeCoords[elementNodes[iElement, :], 1]
        yCoords = nodeCoords[elementNodes[iElement, :], 2]
        plot!(plot1,xCoords, yCoords, color=lc, linestyle=ls, label=false)
    end
end
  
export formStiffness2Dtruss, drawingMesh, stresses2Dtruss, solution
end