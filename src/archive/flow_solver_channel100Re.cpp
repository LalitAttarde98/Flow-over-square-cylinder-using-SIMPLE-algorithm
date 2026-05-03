#include "../include/SIMPLE_solver/SIMPLE_solver.hpp"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm> 
#include <stdexcept>

SimpleSolver::SimpleSolver(int x_grid, int y_grid, double resolution,
                             double U_inlet, double mu_visc, double rho, 
                            double alpha_u, double alpha_v, double alpha_p, double SOR_p)
    : grid_x(x_grid), grid_y(y_grid), resolution(resolution), init_U(U_inlet), mu(mu_visc), rho(rho),
      alpha_u(alpha_u), alpha_v(alpha_v), alpha_p(alpha_p), SOR_p(SOR_p) {
    
    dx = resolution;
    dy = resolution;

    itr = 0;
    res = 10;
    in_flow = init_U * rho * dy * grid_y; // Initial mass flow at inlet

    int nx_tot = grid_x + 2;
    int ny_tot = grid_y + 2;

    auto resize_field = [&](std::vector<std::vector<double>>& field) {
        field.assign(nx_tot, std::vector<double>(ny_tot, 0.0));
    };

    resize_field(p);
    resize_field(p_prime);
    resize_field(u);
    resize_field(v);
    resize_field(u_star);
    resize_field(v_star);
    resize_field(p_star); 
    resize_field(ap_u);
    resize_field(ap_v);
    resize_field(ap_p);
    resize_field(mass_res);
    resize_field(ae); resize_field(aw); resize_field(an); resize_field(as);
    resize_field(u_intp); resize_field(v_intp); resize_field(p_intp); 

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

void SimpleSolver::run(double tolerance) {
    std::cout << std::setprecision(6) << std::fixed << std::endl;
    
    while (res > tolerance) { 
        u_res = solve_u_momentum();
        v_res = solve_v_momentum();
        solve_pressure_correction();
        correct_velocity_and_pressure();
        res = calculate_residual();

        res = res + u_res + v_res;

        if (itr % 10 == 0) {
             std::cout << "Iter: " << itr << " Mass res:" << res << " U_res:" << u_res << " V_res:" << v_res << std::endl;
        }
        itr++;

        if (itr > 10000) {
            std::cout << "Error: Solver did not converge." << std::endl;
            break;
        }
    }

    post_process();
    write_results("results.csv");
}

double SimpleSolver::solve_u_momentum() {
    double fn, fs, fe, fw, de, dw, dn, ds;

    for (int i = 2; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {
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
            an[i][j] = std::max(-fn, std::max(dn - (fn/2), 0.00));
            as[i][j] = std::max(fs, std::max(ds + (fs/2), 0.00));
            ae[i][j] = std::max(-fe, std::max(de - (fe/2), 0.00));
            aw[i][j] = std::max(fw, std::max(dw + (fw/2), 0.00));
    
            ap_u[i][j] = an[i][j] + as[i][j] + ae[i][j] + aw[i][j] + \
                         fn - fs + fe - fw;
            
        }
    }
    // for (int i = 0; i < (grid_x + 2); i++) {
    //     for (int j = 0; j < (grid_y + 2); j++) {
    //     std::cout << u[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // Gauss Seidel
    for (int k = 0; k < 1000; k++) {
        for (int i = 2; i < (grid_x + 1); i++) {
            for (int j = 1; j < (grid_y + 1); j++) {
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
    // for (int i = 0; i < (grid_x + 2); i++) {
    //     for (int j = 0; j < (grid_y + 2); j++) {
    //     std::cout << u[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    float out_flow = 1e-6;
    for (int j = 1; j < grid_y + 1; j++){
        out_flow += (rho * dy * u[grid_x+1][j]);
        }
    //std::cout<<"out_flow: "<< out_flow << " " << (in_flow / out_flow) << std::endl;
    for (int j = 0; j < grid_y + 1; j++){
        u[grid_x+1][j] = (in_flow / out_flow) * u[grid_x+1][j];
        }
    float u_res = 0.0;
    for (int i = 2; i < (grid_x + 1); i++) {
            for (int j = 1; j < (grid_y + 1); j++) {
            u_res += pow(u[i][j] - u_star[i][j], 2.0);
            }
        }
        u_res = sqrt(u_res);

    // for (int i = 0; i < (grid_x + 2); i++) {
    //     for (int j = 0; j < (grid_y + 2); j++) {
    //     std::cout << u[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "U momentum residual: " << u_res << std::endl;
    return u_res;
}

double SimpleSolver::solve_v_momentum() {
    double fn, fs, fe, fw, de, dw, dn, ds;

    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 2; j < grid_y + 1; j++) {
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
            
            an[i][j] = std::max(-fn, std::max(dn - fn / 2, 0.00));
            as[i][j] = std::max(fs, std::max(ds + fs / 2, 0.00));
            ae[i][j] = std::max(-fe, std::max(de - fe / 2, 0.00));
            aw[i][j] = std::max(fw, std::max(dw + fw / 2, 0.00));

            ap_v[i][j] = an[i][j] + as[i][j] + ae[i][j] + aw[i][j] + \
                         fn - fs + fe - fw;

        }
    }

    // Gauss Seidel
    for (int k = 0; k < 1000; k++) {
        for (int i = 1; i < grid_x + 1; i++) {
            for (int j = 2; j < grid_y + 1; j++) {
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
    // for (int i = 0; i < (grid_x + 2); i++) {
    //     for (int j = 0; j < (grid_y + 2); j++) {
    //     std::cout << ap_v[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    float v_res = 0.0;
    for (int i = 1; i < (grid_x + 1); i++) {
        for (int j = 2; j < (grid_y + 1); j++) {
        v_res += pow(v[i][j] - v_star[i][j], 2.0);
        }
    }
    v_res = sqrt(v_res);
    //std::cout << "V momentum residual: " << v_res << std::endl;
    return v_res;
}

void SimpleSolver::solve_pressure_correction() {
    // Discretised coefficients
    for (int i = 1; i < grid_x+1; i++) {
        for (int j = 1; j < grid_y+1; j++) {
            an[i][j] = (rho * dx * dx) / ap_v[i][j + 1];
            as[i][j] = (rho * dx * dx) / ap_v[i][j]; 
            ae[i][j] = (rho * dy * dy) / ap_u[i + 1][j];
            aw[i][j] = (rho * dy * dy) / ap_u[i][j];

            if(i == grid_x) ae[i][j] = 0.0;
            if(i == 1)      aw[i][j] = 0.0;
            if(j == grid_y) an[i][j] = 0.0;
            if(j == 1)      as[i][j] = 0.0;

            ap_p[i][j] = an[i][j] + as[i][j] + ae[i][j] + aw[i][j];
        }
    }
    float residual = 0.0;
    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {
            p_prime[i][j] = 0;
            mass_res[i][j] = (rho*dy*(u[i][j]-u[i+1][j]) + \
                            rho*dx*(v[i][j]-v[i][j+1]));
            residual += pow(mass_res[i][j], 2.0);
        }
    }
    residual = sqrt(residual);
    //std::cout << "Pressure correction residual: " << residual << std::endl;

    ap_p[1][1] = 1.0;
    for (int k = 0; k < 1000; k++) {
        for (int i = 1; i < grid_x+1; i++) {
            for (int j = 1; j < grid_y+1; j++) {
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
            p[i][j] = p[i][j] + alpha_p * p_prime[i][j];
            if (j > 1) {
                v[i][j] = v[i][j] + (dx / ap_v[i][j]) * (p_prime[i][j - 1] - p_prime[i][j]);
                v_star[i][j] = v[i][j];
                }
            if (i > 1) {
                u[i][j] = u[i][j] + (dy / ap_u[i][j]) * (p_prime[i - 1][j] - p_prime[i][j]);
                u_star[i][j] = u[i][j];
            }
        }
    }
}

double SimpleSolver::calculate_residual() {
    double residual = 0.0;
    double mass_balance = 0.0;
    for (int i = 1; i < grid_x + 1; i++) {
        for (int j = 1; j < grid_y + 1; j++) {
            mass_balance = (rho*dy*(u[i][j]-u[i+1][j]) + \
                            rho*dx*(v[i][j]-v[i][j+1]));
            residual += pow(mass_balance, 2.0);
        }
    }

    // for (int i = 0; i < (grid_x + 2); i++) {
    //      for (int j = 0; j < (grid_y + 2); j++) {
    //      std::cout << u[i][j] << " ";
    //      }
    //      std::cout << std::endl;
    //  }
    return std::sqrt(residual);
}


void SimpleSolver::post_process() {
    for (int j = 1; j <= grid_y; j++) {
        for (int i = 1; i <= grid_x; i++) {
            u_intp[i][j] = 0.5 * (u[i][j] + u[i+1][j]);
            v_intp[i][j] = 0.5 * (v[i][j] + v[i][j+1]);
            p_intp[i][j] = p[i][j];
        }
    }
}

void SimpleSolver::write_results(const std::string& filename) const {
    std::ofstream output(filename);
    if (!output.is_open()) {
        std::cerr << "Error: Could not open " << filename << std::endl;
        return;
    }

    output << "x,y,p,u,v" << std::endl;
    for (int j = 1; j <= grid_y; j++) {
        for (int i = 1; i <= grid_x; i++) {
            double x_coord = (i-1) * dx + dx/2.0;
            double y_coord = (j-1) * dy + dy/2.0;
            
            output << x_coord << "," << y_coord << ","
                   << p_intp[i][j] << "," 
                   << u_intp[i][j] << "," 
                   << v_intp[i][j] << std::endl;
        }
    }
    std::cout << "Written to " << filename << std::endl;
}