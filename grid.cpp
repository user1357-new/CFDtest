#include "grid.h"
#include <fstream>
#include <iomanip>
#include <iostream>

ConsVar Grid::primToCons(const State& prim) const {
    double rhoE = prim.p / (config.gamma - 1.0) + 0.5 * prim.rho * prim.u * prim.u;
    return ConsVar(prim.rho, prim.rho * prim.u, rhoE);
}

State Grid::consToPrim(const ConsVar& cons) const {
    double rho = cons.rho;
    double u = cons.rhou / rho;
    double E = cons.rhoE / rho;
    double p = (config.gamma - 1.0) * (E * rho - 0.5 * rho * u * u);
    return State(rho, u, p);
}

Grid::Grid(int n, const Config& cfg) 
    : nx(n), 
      conservativeVars(n), 
      conservativeVarsNew(n),
      config(cfg) {
    initialize();
}

void Grid::initialize() {
    dx = (config.xMax - config.xMin) / nx;
    xCoords.resize(nx);
    for (int i = 0; i < nx; ++i) {
        xCoords[i] = config.xMin + (i + 0.5) * dx;
    }
    
    for (int i = 0; i < nx; ++i) {
        if (xCoords[i] <= 0.0) {
            conservativeVars[i] = primToCons(State(1.0, 0.75, 1.0));
        } else {
            conservativeVars[i] = primToCons(State(0.125, 0.0, 0.1));
        }
    }
    applyBoundaryConditions();
}

void Grid::applyBoundaryConditions() {
    conservativeVars[0] = primToCons(State(1.0, 0.75, 1.0));
    conservativeVars[1] = primToCons(State(1.0, 0.75, 1.0));
    conservativeVars[nx-2] = primToCons(State(0.125, 0.0, 0.1));
    conservativeVars[nx-1] = primToCons(State(0.125, 0.0, 0.1));
}

int Grid::getNx() const { return nx; }
double Grid::getDx() const { return dx; }

std::vector<State> Grid::getPrimitives() const {
    std::vector<State> prim(nx);
    for (int i = 0; i < nx; ++i) {
        prim[i] = consToPrim(conservativeVars[i]);
    }
    return prim;
}

const std::vector<double>& Grid::getXCoords() const { 
    return xCoords; 
}

const std::vector<ConsVar>& Grid::getConservativeVars() const {
    return conservativeVars;
}

void Grid::setConservativeVars(const std::vector<ConsVar>& vars) {
    conservativeVars = vars;
}

void Grid::updateSolution(const std::vector<ConsVar>& flux, double dt) {
    double dt_dx = dt / dx;
    double coeff2 = config.mu * dt / (dx * dx);
    double coeff4 = config.mu4 * dt / (dx * dx);

    for (int i = 1; i < nx-1; ++i) {
        State curr = consToPrim(conservativeVars[i]);
        State left = consToPrim(conservativeVars[i-1]);
        State right = consToPrim(conservativeVars[i+1]);
        ConsVar flux_term = conservativeVars[i] - (flux[i] - flux[i-1]) * dt_dx;
        ConsVar diff2 = conservativeVars[i+1] - conservativeVars[i]*2.0 + conservativeVars[i-1];
        
         double pressure_sensor = fabs(right.p - 2.0 * curr.p + left.p) / (fabs(right.p) + fabs(curr.p) + fabs(left.p) + 1e-10);
        ConsVar visc2_term, visc4_term;
        visc4_term.rho = 0.0; visc4_term.rhou = 0.0; visc4_term.rhoE = 0.0;
        visc2_term = diff2 * coeff2;
         if (pressure_sensor > 0.05) { 
            // --- 激波区 ---
            visc2_term = diff2 * (coeff2 * 2.0);
        } else {
            // --- 光滑区 ---
            //visc2_term = diff2 * coeff2;
        if (i >= 2 && i <= nx - 3 ){
        ConsVar diff4 = conservativeVars[i+2] 
                      - conservativeVars[i+1]*4.0 
                      + conservativeVars[i]*6.0 
                      - conservativeVars[i-1]*4.0 
                      + conservativeVars[i-2];
        visc4_term =diff4 * coeff4;
    }
}
         conservativeVarsNew[i] = flux_term + visc2_term + visc4_term;
    
    }
    
    for (int i = 1; i < nx - 1; ++i) {
        if (conservativeVarsNew[i].rho < 1e-6) {
            conservativeVarsNew[i].rho = 1e-6;
        }
    }
    std::swap(conservativeVars, conservativeVarsNew);
    applyBoundaryConditions();
}

void Grid::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "无法打开文件: " << filename << std::endl;
        return;
    }

    auto prim = getPrimitives();
    file << std::fixed << std::setprecision(6);
    for (int i = 0; i < nx; ++i) {
        file << xCoords[i] << " "
             << prim[i].rho << " "
             << prim[i].u << " "
             << prim[i].p << std::endl;
    }

    file.close();
    std::cout << "save " << filename << std::endl;
}