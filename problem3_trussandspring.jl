##################################################################
###########        FEM Lab Course -by Saurabh Singh    ###########
##################################################################
push!(LOAD_PATH,"./")

using LinearAlgebra
using Plots
# using functions

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
function drawingMesh(plot1, nodeCoords, elementNodes,N_elems, lc,ls)
  for iElement in 1:N_elems
      xCoords = nodeCoords[elementNodes[iElement, :], 1]
      yCoords = nodeCoords[elementNodes[iElement, :], 2]
      plot!(plot1,xCoords, yCoords, color=lc, linestyle=ls, label=false)
  end
end

function main()
  # Write a program to solve a general 2D truss problem using FEM
  # Daryl Logan Example 3.5
  # antonio ferreira 2008
  # Modified by Ahmed Rashed
  # Units used are N-mm

  # Material properties
  E_vec = 210000*ones(2) #11 element vector with value 70e9
  A_vec = 500*ones(2)
  k     = 2000 #N/mm

  # Generation of coordinates and connectivities
  nodesCoords = [
    0. 0.;
    -5000.0*cos(pi/4) 5000.0*sin(pi/4);
    -10000.0 0.
  ]
  elementNodes = [
    1 2;
    1 3;
    1 4
  ]

  numberElements = size(elementNodes, 1)
  numberNodes = size(nodesCoords, 1) + 1 #+1 for spring node

  xx = nodesCoords[:,1]
  yy = nodesCoords[:,2]

  # GDof: total number of degrees of freedom
  GDof = 2 * (numberNodes)

  # Displacement vector
  D_col = fill(NaN, GDof)
  # Force vector
  F_col = fill(NaN, GDof)

  #Applied Force 
  F_col[2] = -25000
  F_col[1] = 0
  # F_col[8] = 0

  # Boundary conditions
  prescribedDof = 3:8
  for idof in prescribedDof 
    D_col[idof] = 0.0 
  end
  
  # K_assembly = zeros(GDof,GDof)
  # Assembly stiffness matrix
  K_assembly = formStiffness2Dtruss(GDof, numberElements-1, elementNodes, nodesCoords, E_vec, A_vec)

  # display(K_assembly)
  #Spring stiffness
  K_assembly[[2,7], [2,7]] += k*[1 -1; -1 1]

  # Solution
  D_col, F_col = solution(prescribedDof, K_assembly, D_col, F_col)

  ### POST PROCESSING
  # Deformation components
  us = 1:2:(2 * (numberNodes-1) - 1)
  vs = 2:2:(2 * (numberNodes-1))

  # Extract the x and y components of the displacement
  XX = D_col[us]
  YY = D_col[vs]

  # Calculate the maximum norm of the displacements
  dispNorm = maximum(sqrt.(XX.^2 .+ YY.^2))

  # Scale factor for displacements
  scaleFact =  2*dispNorm

  # Adjust node coordinates with displacements for visualization
  adjustedCoords = nodesCoords .+ scaleFact * hcat(XX, YY)

  # Initial plot setup
  plot1 = plot(aspect_ratio=:equal)

  # Draw the original mesh
  drawingMesh(plot1, nodesCoords, elementNodes,numberElements-1,:black, :dot)

  # Draw the deformed mesh
  drawingMesh(plot1,adjustedCoords, elementNodes,numberElements-1,:red,:solid)

  # Save the plot
  savefig("problem3_truss&spring.png")
  # Stress at elements
  stress = stresses2Dtruss(numberElements-1, elementNodes, nodesCoords, D_col, E_vec)

  # Display Result
  println("\nStiffness Matrix: ")
  display(K_assembly)
  println("\nDisplacement Vector: ")
  display(D_col)
  println("\nForce Vector: ")
  display(F_col)
  println("\nStress: ")
  display(stress)
      
end

main()