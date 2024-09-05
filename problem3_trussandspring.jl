##################################################################
###########        FEM Lab Course -by Saurabh Singh    ###########
##################################################################
push!(LOAD_PATH,"./")

using LinearAlgebra
using Plots
using functions

function main()
  # Write a program to solve a general 2D truss problem using FEM
  # Daryl Logan Example 3.5
  # antonio ferreira 2008
  # Modified by Ahmed Rashed
  # Units used are N-mm

  # Material properties
  E_vec = 70000*ones(11) #11 element vector with value 70e9
  A_vec = 300*ones(11)

  # Generation of coordinates and connectivities
  nodesCoords = [
    0. 0.;
    -5000.0*cos(pi/4) 5000.0*sin(pi/4);
    -10000.0 0.
  ]
  numberNodes = size(nodesCoords, 1) 

  elementNodes = [
    1 2;
    1 3;
    1 4
  ]
  numberElements = size(elementNodes, 1)

  # GDof: total number of degrees of freedom
  GDof = 2 * (numberNodes) #+1 for spring

  # Assembly stiffness matrix
  K_assembly = formStiffness2Dtruss(GDof, numberElements, elementNodes, nodesCoords, E_vec, A_vec)

  #Spring stiffness
  K_assembly[[2,7], [2,7]] += 2000*[1 -1; -1 1]

  # Boundary conditions
  prescribedDof = 3:8

  # Force vector
  F_col = fill(NaN, GDof)
  F_col[2] = -25e3

  # Displacement vector
  D_col = fill(NaN, GDof)

  # Solution
  D_col, F_col = solution(prescribedDof, K_assembly, D_col, F_col)

  ### POST PROCESSING
  # Deformation components
  us = 1:2:(2 * numberNodes - 1)
  vs = 2:2:(2 * numberNodes)

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
  drawingMesh(plot1, nodesCoords, elementNodes,:black, :dot)

  # Draw the deformed mesh
  drawingMesh(plot1,adjustedCoords, elementNodes,:red,:solid)

  # Save the plot
  savefig("problem3_truss&spring.png")
  # Stress at elements
  stress = stresses2Dtruss(numberElements, elementNodes, nodesCoords, D_col, E_vec)

  # Display Result
  println("\nDisplacement Vector: ")
  display(D_col)
  println("\nForce Vector: ")
  display(F_col)
  println("\nStress: ")
  display(stress)
      
end

main()