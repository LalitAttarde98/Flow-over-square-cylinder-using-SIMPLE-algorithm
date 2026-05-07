# Flow-over-square-cylinder-using-SIMPLE-algorithm

The present work models a 2D flow around a square cylinder
using a finite volume approach with Semi Implicit method
for Pressure Linked Equations. The analysis focuses on the 
effect of Reynolds number, as determined by the input velocity
parameters, on the flow around the square cylinder. Velocity profiles 
along with an integral drag coefficient are investigated for
laminar flows, and validated with the literature[3]. Refer to the technical ![report](doc/Project_Report.pdf) for implementation details and validation results.

![Results](doc\results_100Re_flow_fields.png)

# How to run the solver
1) Define flow, grid, solver params, and position, size & number of obstacles in src/main.cpp
2) Using cmake to build project
   `mkdir build && cd build # Create build folder
   `cmake ..                # Configure
   `cmake --build .         # Build
   `./run_solver            # Run

References:
1. Suhas V. Patankar, Textbook, Numerical Heat Transfer and
Flud Flow
2. Fadl Moukalled, Luca Mangani, Marwan Darwish (2016)
Implementation of boundary conditions in the finite-volume
pressure-based method—Part I: Segregated solvers,
Numerical Heat Transfer, Part B: Fundamentals, 69:6, 534-562,
DOI: 10.1080/10407790.2016.1138748
3. Breuer, M, Bernsdorf, J., ans T.Zeiser, F. D., 1999. 
“Accurate computations of the laminar flow past a square cylinder
 based on two different methods: lattice boltzmann and
finite-volume
4. H K Versteeg, W Malalasekara, Textbook, An Introduction
to computational Fluid Dynamics, Finite Volume Method
