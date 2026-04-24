
#include "fvs_solver.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

FVSSolver::FVSSolver(Grid& g, const Config& cfg, int spatialOrder, int timeMethod) 
    : grid(g), config(cfg), order(spatialOrder), timeScheme(timeMethod) {
}

double FVSSolver::computeWaveSpeed(const State& state) const {
    double c = sqrt(config.gamma * state.p / state.rho);
    return fabs(state.u) + c;
}

ConsVar FVSSolver::computeStegerWarmingFlux(const State& state) const {
    double c = sqrt(config.gamma * state.p / state.rho);
    double rho = state.rho;
    double u = state.u;
    double p = state.p;
    
    // 守恒变量
    double rhoE = p / (config.gamma - 1.0) + 0.5 * rho * u * u;
    double H = (config.gamma / (config.gamma - 1.0)) * p / rho + 0.5 * u * u;
    
    // Steger-Warming 特征值分裂
    double a1 = std::max(0.0, u);
    double a2 = std::max(0.0, u + c);
    double a3 = std::max(0.0, u - c);
    
    // 正通量部分
    double f1 = (a1 + (a2 + a3) / (2.0 * config.gamma)) * rho;
    double f2 = (a1 + (a2 + a3) / (2.0 * config.gamma)) * rho * u 
                + (a2 - a3) * rho * c / (2.0 * config.gamma);
    double f3 = (a1 + (a2 + a3) / (2.0 * config.gamma)) * rhoE 
                + (a2 - a3) * rho * c * H / (2.0 * config.gamma);
    
    return ConsVar(f1, f2, f3);
}

ConsVar FVSSolver::computeRusanovFlux(const State& left, const State& right) const {
    ConsVar leftCons = grid.primToCons(left);
    ConsVar rightCons = grid.primToCons(right);
    
    // 物理通量
    auto physFlux = [](const State& s, double gamma) -> ConsVar {
        double E = s.p / ((gamma - 1.0) * s.rho) + 0.5 * s.u * s.u;
        return ConsVar(
            s.rho * s.u,
            s.rho * s.u * s.u + s.p,
            s.u * (s.rho * E + s.p)
        );
    };
    
    ConsVar leftFlux = physFlux(left, config.gamma);
    ConsVar rightFlux = physFlux(right, config.gamma);
    
    // 最大波速
    double alpha = std::max(computeWaveSpeed(left), computeWaveSpeed(right));
    
    // Rusanov 通量
    ConsVar flux;
    flux.rho  = 0.5 * (leftFlux.rho + rightFlux.rho) - 0.5 * alpha * (rightCons.rho - leftCons.rho);
    flux.rhou = 0.5 * (leftFlux.rhou + rightFlux.rhou) - 0.5 * alpha * (rightCons.rhou - leftCons.rhou);
    flux.rhoE = 0.5 * (leftFlux.rhoE + rightFlux.rhoE) - 0.5 * alpha * (rightCons.rhoE - leftCons.rhoE);
    
    return flux;
}

double FVSSolver::minmod(double a, double b) const {
    if (a * b > 0.0) {
        return (fabs(a) < fabs(b)) ? a : b;
    }
    return 0.0;
}

State FVSSolver::minmod(const State& a, const State& b) const {
    return State(
        minmod(a.rho, b.rho),
        minmod(a.u, b.u),
        minmod(a.p, b.p)
    );
}

State FVSSolver::reconstructState(int i, const std::vector<State>& primitives, bool leftSide) const {
    int nx = grid.getNx();
    
    if (order == 1) {
        return primitives[i];
    }
    
    const State& ui = primitives[i];
    const State& ui_minus1 = primitives[(i-1 >= 0) ? i-1 : i];
    const State& ui_plus1 = primitives[(i+1 < nx) ? i+1 : i];

    State dL = ui - ui_minus1;
    State dR = ui_plus1 - ui;
    
    State reconstructed = ui + minmod(dL, dR) * 0.5;

    return reconstructed;
}

double FVSSolver::computeTimeStep() {
    auto primitives = grid.getPrimitives();
    double maxSpeed = 0.0;
    double dx = grid.getDx();
    
    for (const auto& state : primitives) {
        maxSpeed = std::max(maxSpeed, computeWaveSpeed(state));
    }
    
    return config.cfl * dx / maxSpeed;
}

void FVSSolver::eulerStep(std::vector<ConsVar>& U_new, const std::vector<ConsVar>& U_old, double dt) {
    int nx = grid.getNx();
    
    grid.setConservativeVars(U_old);
    auto primitives = grid.getPrimitives();
    
    std::vector<ConsVar> fluxes(nx);
    
    for (int i = 1; i < nx; ++i) {
        State leftState, rightState;
        
        if (order == 1) {
            leftState = primitives[i-1];
            rightState = primitives[i];
        } else {
            leftState = reconstructState(i-1, primitives, false);
            rightState = reconstructState(i, primitives, true);
        }
        
        fluxes[i] = computeRusanovFlux(leftState, rightState);
    }
    
    double dt_dx = dt / grid.getDx();
    for (int i = 1; i < nx - 1; ++i) {
        U_new[i] = U_old[i] - (fluxes[i+1] - fluxes[i]) * dt_dx;
    }
    
    U_new[0] = U_new[1];
    U_new[nx-1] = U_new[nx-2];
}

void FVSSolver::rk3Step(std::vector<ConsVar>& U_final, const std::vector<ConsVar>& U_old, double dt) {
    int nx = grid.getNx();
    
    // TVD-RK3 第一阶段
    std::vector<ConsVar> U1(nx);
    eulerStep(U1, U_old, dt);
    
    // TVD-RK3 第二阶段
    std::vector<ConsVar> U1_tilde(nx);
    eulerStep(U1_tilde, U1, dt);
    std::vector<ConsVar> U2(nx);
    for (int i = 0; i < nx; ++i) {
        U2[i] = U_old[i] * 0.75 + U1_tilde[i] * 0.25;
    }
    
    // TVD-RK3 第三阶段
    std::vector<ConsVar> U2_tilde(nx);
    eulerStep(U2_tilde, U2, dt);
    for (int i = 0; i < nx; ++i) {
        U_final[i] = U_old[i] / 3.0 + U2_tilde[i] * 2.0 / 3.0;
    }
}

void FVSSolver::evolve() {
    double currentTime = 0.0;
    int step = 0;
    int nx = grid.getNx();
    
    std::cout << "=== FVS Solver (Space=" << order 
              << ", Time=" << (timeScheme == 1 ? "Euler" : "RK3") << ") ===" << std::endl;
    
    while (currentTime < config.tEnd) {
        double dt = computeTimeStep();
        if (currentTime + dt > config.tEnd) {
            dt = config.tEnd - currentTime;
        }
        
        auto U = grid.getConservativeVars();
        std::vector<ConsVar> U_new(nx);
        
        if (timeScheme == 1) {
            eulerStep(U_new, U, dt);
        } else {
            rk3Step(U_new, U, dt);
        }
        
        grid.setConservativeVars(U_new);
        
        currentTime += dt;
        step++;
        
        if (step % 10 == 0 || currentTime >= config.tEnd) {
            std::cout << "Step: " << step 
                     << " | Time: " << std::fixed << std::setprecision(4) << currentTime
                     << " | dt: " << dt << std::endl;
        }
    }
    
    std::cout << "FVS 推进完成！总步数：" << step << std::endl;
}