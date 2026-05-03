#include "../include/SIMPLE_solver/SIMPLE_solver.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm> 
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265359
#endif

SimpleSolver::SimpleSolver(int x_grid, int y_grid, float resolution, int geom_height,
                             float U_inlet, float mu_visc, float rho, 
                            float alpha_u, float alpha_v, float alpha_p, float SOR_p)
    : grid_x(x_grid), grid_y(y_grid), resolution(resolution), geom_height(geom_height),
      init_U(U_inlet), mu(mu_visc), rho(rho),
      alpha_u(alpha_u), alpha_v(alpha_v), alpha_p(alpha_p), SOR_p(SOR_p) {
    
    dx = resolution;
    dy = resolution;

    in_flow = init_U * rho * dy * grid_y; 

    // Include ghost cells (0 and grid+1)
    nx_tot = grid_x + 2;
    ny_tot = grid_y + 2;

    auto resize_field = [&](std::vector<std::vector<float>>& field) {
        field.assign(nx_tot, std::vector<float>(ny_tot, 0.0f));
    };

    resize_field(p);
    resize_field(p_prime); 
    resize_field(u);
    resize_field(v);
    resize_field(u_star);
    resize_field(v_star);
    resize_field(ap_u);
    resize_field(ap_v);
    resize_field(ap_p);
    resize_field(mass_res);
    resize_field(ae); resize_field(aw); resize_field(an); resize_field(as); 

    initialize_fields();
}

void SimpleSolver::initialize_fields() {

    for (int i = 0; i < grid_x + 2; i++) {
        for (int j = 0; j < grid_y + 2; j++) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
            p[i][j] = 0.0;
            u_star[i][j] = 0.0;
            v_star[i][j] = 0.0;
            p_prime[i][j] = 0.0;

            if (i == 1 || i == grid_x + 1) {
                u[i][j] = init_U; 
                u_star[i][j] = init_U;
            }
        }
    }
}

void SimpleSolver::run(float tolerance, int total_interval, float alpha_u_min, float alpha_v_min, float alpha_p_min) {
    std::cout << std::setprecision(6) << std::fixed << std::endl;
    float alpha_p_max = alpha_p;
    float alpha_u_max = alpha_u;
    float alpha_v_max = alpha_v;
    float total_error = 1e6;
    int step = 0;
    while (total_error > tolerance) { 
        float u_res = solve_u_momentum();
        float v_res = solve_v_momentum();
        solve_pressure_correction();
        correct_velocity_and_pressure();
        float m_res = calculate_residual();
        total_error = m_res + u_res + v_res;
        // Consine annealing for relaxation factors
        alpha_p = alpha_p_min + (alpha_p_max - alpha_p_min) *  \
            (0.5 * (1.0 + std::cos(M_PI * step / total_interval)));
        alpha_u = alpha_u_min + (alpha_u_max - alpha_u_min) *  \
            (0.5 * (1.0 + std::cos(M_PI * step / total_interval)));
        alpha_v = alpha_v_min + (alpha_v_max - alpha_v_min) *  \
            (0.5 * (1.0 + std::cos(M_PI * step / total_interval)));
        // Stepwise reduction of relaxation factors
        if (step == 50) alpha_u_max = alpha_u_max / 2.0;
        if (step == 50) alpha_v_max = alpha_v_max / 2.0; 

        if (step % 10 == 0) {
             std::cout << "Iter: " << step << " Mass res:" << m_res << " U_res:" << u_res << " V_res:" << v_res << \
              " alpha_p:" << alpha_p << " alpha_u:" << alpha_u << " alpha_v:" << alpha_v << " Total error:" << total_error << std::endl;
        }
        step++;

        if (step > total_interval) {
            std::cout << "Error: Solver did not converge." << std::endl;
            break;
        }
    }
    calculate_drag_coefficient();
}

void SimpleSolver::create_geometry(const std::vector<int>& x_center, const std::vector<int>& y_center, int geom_width, int geom_height) {

    geometry.assign(grid_x + 2, std::vector<bool>(grid_y + 2, false));
    size_t num_obstacles = x_center.size();

    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {
            bool is_obstacle = false;
            for (int k = 0; k < num_obstacles; k++) {
                if (i >= x_center[k] - geom_width/2 && i <= x_center[k] + geom_width/2 && \
                    j >= y_center[k] - geom_height/2 && j <= y_center[k] + geom_height/2) {
                    is_obstacle = true;
                    break; 
                }
            }
            geometry[i][j] = is_obstacle;
            if (is_obstacle) {
                u[i][j] = 0.0;
                v[i][j] = 0.0;
                u_star[i][j] = 0.0;
                v_star[i][j] = 0.0;
                p_prime[i][j] = 0.0;
                p[i][j] = 0.0;
            }
        }
    }
}

float SimpleSolver::solve_u_momentum() {
    float fn, fs, fe, fw, de, dw, dn, ds;

    for (int i = 2; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {

            if (geometry[i][j] || geometry[i-1][j]) continue;
            // Convective fluxes
            fn = 0.5 * dx * rho * (v_star[i - 1][j + 1] + v_star[i][j + 1]);
            fs = 0.5 * dx * rho * (v_star[i - 1][j] + v_star[i][j]);
            fe = 0.5 * dy * rho * (u_star[i][j] + u_star[i + 1][j]);
            fw = 0.5 * dy * rho * (u_star[i - 1][j] + u_star[i][j]);

            de = (mu * dy) / dx;
            dw = (mu * dy) / dx;
            dn = (mu * dx) / dy;
            ds = (mu * dx) / dy;

            if (j == grid_y) dn = 2.0 * dn;  // Half cell at north
            if (j == 1) ds = 2.0 * ds; // Half cell at south
            
            // Hybrid Scheme
            an[i][j] = std::max(-fn, std::max(dn - (fn/2), 0.0f));
            as[i][j] = std::max(fs, std::max(ds + (fs/2), 0.0f));
            ae[i][j] = std::max(-fe, std::max(de - (fe/2), 0.0f));
            aw[i][j] = std::max(fw, std::max(dw + (fw/2), 0.0f));
    
            ap_u[i][j] = an[i][j] + as[i][j] + ae[i][j] + aw[i][j] + \
                         fn - fs + fe - fw;
            
        }
    }
    // Gauss seidel solver
    for (int k = 0; k < 500; k++) {
        for (int i = 2; i < (grid_x + 1); i++) {
            for (int j = 1; j < (grid_y + 1); j++) {
    
                if (geometry[i][j] || geometry[i-1][j]) {
                    u[i][j] = 0.0; 
                    continue;
                }

                u[i][j] = (1.0 - alpha_u) * u_star[i][j] + \
                        (alpha_u / ap_u[i][j]) * 
                        ((an[i][j] * u[i][j + 1]) + \
                        (as[i][j] * u[i][j - 1]) + \
                        (ae[i][j] * u[i + 1][j]) + \ 
                        (aw[i][j] * u[i - 1][j]) + \
                        (dy * (p[i - 1][j] - p[i][j])));
        
                if (std::isnan(u[i][j])) {
                        std::cout << i << j << "t";
                        throw std::runtime_error("NaN detected in u-momentum!");
                }
            }
        }
        for (int j = 0; j < grid_y + 2; j++) {
            u[grid_x+1][j] = u[grid_x][j]; // zero gradient outlet
        }
    }

    float out_flow = 1e-6;
    for (int j = 1; j < grid_y + 1; j++){
        out_flow += (rho * dy * u[grid_x+1][j]);
        }
    //std::cout<<"out_flow: "<< out_flow << " " << (in_flow / out_flow) << std::endl;
    if (out_flow > 1e-3){
        for (int j = 0; j < grid_y + 1; j++){
        u[grid_x+1][j] = (in_flow / out_flow) * u[grid_x+1][j]; // mass flow correction at outlet
        }
    }
    float u_res = 0.0;
    for (int i = 2; i < (grid_x + 1); i++) {
            for (int j = 1; j < (grid_y + 1); j++) {
            u_res += pow(u[i][j] - u_star[i][j], 2.0);
            }
        }
        u_res = sqrt(u_res);

    return u_res;
}

float SimpleSolver::solve_v_momentum() {
    float fn, fs, fe, fw, de, dw, dn, ds;

    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 2; j < grid_y + 1; j++) {
            if (geometry[i][j] || geometry[i][j-1]) continue; 

            fn = 0.5 * dx * rho * (v_star[i][j] + v_star[i][j + 1]);
            fs = 0.5 * dx * rho * (v_star[i][j - 1] + v_star[i][j]);
            fe = 0.5 * dy * rho * (u_star[i + 1][j - 1] + u_star[i + 1][j]);
            fw = 0.5 * dy * rho * (u_star[i][j - 1] + u_star[i][j]);

            de = (mu * dy) / dx;
            dw = (mu * dy) / dx;
            dn = (mu * dx) / dy;
            ds = (mu * dx) / dy;

            if (i == grid_x) de = 2.0 * de;  // Half cell at east
            if (i == 1) dw = 2.0 * dw; // Half cell at west
            
            an[i][j] = std::max(-fn, std::max(dn - fn / 2, 0.0f));
            as[i][j] = std::max(fs, std::max(ds + fs / 2, 0.0f));
            ae[i][j] = std::max(-fe, std::max(de - fe / 2, 0.0f));
            aw[i][j] = std::max(fw, std::max(dw + fw / 2, 0.0f));

            ap_v[i][j] = an[i][j] + as[i][j] + ae[i][j] + aw[i][j] + \
                         fn - fs + fe - fw;

        }
    }

    for (int k = 0; k < 500; k++) {
        for (int i = 1; i < grid_x + 1; i++) {
            for (int j = 2; j < grid_y + 1; j++) {
                if (geometry[i][j] || geometry[i][j-1]) {
                    v[i][j] = 0.0; 
                    continue;
                }

                v[i][j] = ((1 - alpha_v) * v_star[i][j]) + \
                        (alpha_v / ap_v[i][j]) * \
                        ((an[i][j] * v[i][j + 1]) + \
                        (as[i][j] * v[i][j - 1]) + \
                        (ae[i][j] * v[i + 1][j]) + \
                        (aw[i][j] * v[i - 1][j]) + \
                        (dx * (p[i][j-1] - p[i][j])));
                if (std::isnan(v[i][j])) {
                    std::cout << i << " " << j << "\t";
                    throw std::runtime_error("NaN detected in v-momentum!");
                }
            }
        }
    }
    float v_res = 0.0;
    for (int i = 1; i < (grid_x + 1); i++) {
        for (int j = 2; j < (grid_y + 1); j++) {
        v_res += pow(v[i][j] - v_star[i][j], 2.0);
        }
    }
    v_res = sqrt(v_res);

    return v_res;
}

void SimpleSolver::solve_pressure_correction() {

    for (int i = 1; i < grid_x+1; i++) {
        for (int j = 1; j < grid_y+1; j++) {
            an[i][j] = (rho * dx * dx) / ap_v[i][j + 1];
            as[i][j] = (rho * dx * dx) / ap_v[i][j]; 
            ae[i][j] = (rho * dy * dy) / ap_u[i + 1][j];
            aw[i][j] = (rho * dy * dy) / ap_u[i][j];

            if (i == grid_x) ae[i][j] = 0.0;
            if (i == 1)      aw[i][j] = 0.0;
            if (j == grid_y) an[i][j] = 0.0;
            if (j == 1)      as[i][j] = 0.0;

            if (geometry[i][j] || geometry[i][j+1]) an[i][j] = 0.0;
            if (geometry[i][j] || geometry[i][j-1]) as[i][j] = 0.0;                   
            if (geometry[i][j] || geometry[i+1][j]) ae[i][j] = 0.0;
            if (geometry[i][j] || geometry[i-1][j]) aw[i][j] = 0.0;

            ap_p[i][j] = an[i][j] + as[i][j] + ae[i][j] + aw[i][j];

            p_prime[i][j] = 0;
            mass_res[i][j] = (rho*dy*(u[i][j]-u[i+1][j]) + \
                             rho*dx*(v[i][j]-v[i][j+1]));
        }
    }
    
    //ap_p[1][1] = 1.0; // to avoid singular matrix
    for (int k = 0; k < 500; k++) {
        for (int i = 1; i < grid_x+1; i++) {
            for (int j = 1; j < grid_y+1; j++) {
                if (geometry[i][j]) {
                    p_prime[i][j] = 0.0; 
                    continue;
                }
                p_prime[i][j] = p_prime[i][j] + (SOR_p/ap_p[i][j]) * \
                        (an[i][j] * p_prime[i][j+1] + \
                        as[i][j] * p_prime[i][j-1] + \
                        ae[i][j] * p_prime[i+1][j] + \
                        aw[i][j] * p_prime[i-1][j] + \
                        mass_res[i][j] - (ap_p[i][j] * p_prime[i][j])); 
            }
        }
    }
}

void SimpleSolver::correct_velocity_and_pressure() {

    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {
            if (geometry[i][j]) continue;

            p[i][j] = p[i][j] + alpha_p * p_prime[i][j];
            if ((j > 1) &&  !(geometry[i][j] || geometry[i][j-1])) {
                v[i][j] = v[i][j] + (dx / ap_v[i][j]) * (p_prime[i][j - 1] - p_prime[i][j]);
                v_star[i][j] = v[i][j];
                }
            if ((i > 1) && !(geometry[i][j] || geometry[i-1][j])) {
                u[i][j] = u[i][j] + (dy / ap_u[i][j]) * (p_prime[i - 1][j] - p_prime[i][j]);
                u_star[i][j] = u[i][j];
            }
        }
    }
}

float SimpleSolver::calculate_residual() {
    float residual = 0.0;
    float mass_balance = 0.0;
    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {
            mass_balance = (rho*dy*(u[i][j]-u[i+1][j]) + \
                            rho*dx*(v[i][j]-v[i][j+1]));
            residual += pow(mass_balance, 2.0);
        }
    }
    float m_res = std::sqrt(residual);
    return m_res;
}

void SimpleSolver::calculate_drag_coefficient() {
    float P_front = 0.0;
    float P_back = 0.0;
    float F_viscous = 0.0;
    float area = 0.0;
    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {

            if (geometry[i+1][j] || geometry[i][j]) P_back +=  p[i+1][j];
            if (geometry[i-1][j] || geometry[i][j]) P_front +=  p[i-1][j];
            if (geometry[i][j+1] || geometry[i][j]){
                F_viscous += (mu * dx * 0.5 * ((u[i][j] + u[i+1][j]) / dy));
            }
            if (geometry[i][j-1] || geometry[i][j]){
                F_viscous += (mu * dx * 0.5 * ((u[i][j] + u[i+1][j]) / dy));
            }
        }
    }
    float drag_coefficient = ((P_front - P_back) * dy + F_viscous \
            ) / (0.5 * rho * init_U * init_U * float(geom_height) * resolution);

    std::cout << "Drag Coefficient: " << drag_coefficient << ", Pressure force (N): " \
         << (P_front - P_back) * dy << ", Viscous force (N): " << F_viscous << std::endl;
}

void SimpleSolver::write_results(const std::string& filename) const {
    std::ofstream output(filename);
    if (!output.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }

    output << "x,y,p,u,v" << std::endl;
    float x_coord, y_coord, u_intp, v_intp, p_intp;
    for (int j = 1; j <= grid_y; j++) {
        for (int i = 1; i <= grid_x; i++) {
            x_coord = (i-1) * dx + dx/2.0f;
            y_coord = (j-1) * dy + dy/2.0f;
            // Interpolate u and v to cell centers of pressure
            u_intp = 0.5 * (u[i][j] + u[i+1][j]);
            v_intp = 0.5 * (v[i][j] + v[i][j+1]);
            p_intp = p[i][j];

            if (geometry[i][j] || geometry[i-1][j]) u_intp = 0.0;
            if (geometry[i][j] || geometry[i][j-1]) v_intp = 0.0;
            if (geometry[i][j]) p_intp = 0.0;
            
            output << x_coord << "," << y_coord << ","
                   << p_intp << "," 
                   << u_intp << "," 
                   << v_intp << std::endl;
        }
    }
    std::cout << "Written to " << filename << std::endl;
}