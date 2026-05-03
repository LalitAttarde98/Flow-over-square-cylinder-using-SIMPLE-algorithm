#include "../include/SIMPLE_solver/SIMPLE_solver.hpp" 
#include <iostream>
#include <stdexcept> 
#include <vector>

int main() {
    try {
        float resolution = 0.1; // 10 cm grid resolution
        int y_grid = 40;         // No of grid cells in y direction
        int x_grid = y_grid * 4; // No of grid cells in x direction

        // Flow parameters
        float U_inlet = 1.0; // inlet velocity m/sec
        float Re = 100.0;    // Reynolds number
        float rho = 1000.0; // Water density kg/m^3

        // Solver parameters
        float  tolerance = 1e-4;
        float SOR_p = 1.5;
        float alpha_u = 0.8;
        float alpha_u_min = 0.00001;
        float alpha_v = 0.8;
        float alpha_v_min = 0.00001;
        float alpha_p = 0.1;
        float alpha_p_min = 0.00001;
        int total_interval = 1200;

        // Define geometry (x_center, y_center) of obstacle inside confined flow
        // blocking ratio: 1/4
        std::vector<int> x_center = {x_grid / 4,  (5 * x_grid) / 8};
        int geom_width = (1.0/4.0) * (y_grid);
        std::vector<int> y_center = {(3 * y_grid) / 4, y_grid / 4};
        int geom_height = (1.0/4.0) * (y_grid); 

        // Dynamic viscosity 
        float mu_visc = (U_inlet * geom_height * resolution * rho) / Re; 

        std::cout << "SIMPLE algorithm for staggered grid" << std::endl;
        std::cout << "Flow over square cylinder in the duct" << std::endl;
        std::cout << "Grid size " << x_grid << "*" << y_grid << std::endl;
        std::cout << "Inlet velocity = " << U_inlet << " m/s dynamic viscosity = " << mu_visc << " Kg/m.s" << std::endl;
        std::cout << "Reynold Number =  " << Re << std::endl;

        SimpleSolver solver(x_grid, y_grid, resolution, geom_height,
                             U_inlet, mu_visc, rho,
                            alpha_u, alpha_v, alpha_p, SOR_p);

        solver.create_geometry(x_center, y_center, geom_width, geom_height);

        solver.run(tolerance, total_interval, alpha_u_min, alpha_v_min, alpha_p_min);

        solver.write_results("../logs/results_100Re.csv"); // give the path for output file

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1; 
    }

    return 0; 
}