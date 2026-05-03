#pragma once

#include <vector>
#include <string>

using Field = std::vector<std::vector<float>>;

class SimpleSolver {
public:
    SimpleSolver(int x_grid, int y_grid, float resolution, int geom_height,
                             float U_inlet, float mu_visc, float rho,
                            float alpha_u, float alpha_v, float alpha_p, float SOR_p);

    void create_geometry(const std::vector<int>& x_center, const std::vector<int>& y_center, int geom_width, int geom_height);
    void run(float tolerance, int total_interval, float alpha_u_min, float alpha_v_min, float alpha_p_min);
    void write_results(const std::string& filename) const;

private:
    // Physics Constants
    float rho;        
    float mu;
    float init_U;
    float Re;
    
    // Grid & Geometry
    float resolution;
    int grid_x, grid_y;
    int nx_tot, ny_tot;
    int geom_width, geom_height;
    float dx, dy;

    // Relaxation Factors
    float tolerance;
    float SOR_p;
    float alpha_u, alpha_v,alpha_p;
    float alpha_u_min, alpha_v_min, alpha_p_min;
    int total_interval;

    // Solver State
    float m_res, u_res, v_res;
    float in_flow;

    // Primary flow Variables
    Field u;
    Field v;
    Field p;

    // Intermediate flow Variables
    Field p_prime; 
    Field u_star;
    Field v_star;

    // Coefficients for Discretized Equations
    Field ap_u;
    Field ap_v;
    Field ap_p;
    Field an;
    Field as;
    Field ae;
    Field aw;

    Field mass_res;
    std::vector<std::vector<bool>> geometry;

    void initialize_fields();
    float solve_u_momentum();
    float solve_v_momentum();
    void solve_pressure_correction();
    void correct_velocity_and_pressure();
    float calculate_residual();
    void calculate_drag_coefficient();
    
};