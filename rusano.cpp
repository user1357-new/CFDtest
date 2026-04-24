#include "rusano.h"
#include <algorithm>
#include <iostream>
#include <iomanip>

RusanovSolver::RusanovSolver(Grid& g, const Config& cfg) 
    : grid(g), config(cfg) {
}

ConsVar RusanovSolver::computePhysicalFlux(const State& state) const {
    double rho = state.rho;
    double u = state.u;
    double p = state.p;
    double E = p / ((config.gamma - 1.0) * rho) + 0.5 * u * u;
    return ConsVar(rho * u,
                  rho * u * u + p,
                  u * (rho * E + p));
}

double RusanovSolver::computeWaveSpeed(const State& state) const {
    double soundSpeed = sqrt(config.gamma * state.p / state.rho);
    return fabs(state.u) + soundSpeed;
}

ConsVar RusanovSolver::computeRusanovFlux(const State& leftState, const State& rightState) {
    ConsVar leftCons = grid.primToCons(leftState);
    ConsVar rightCons = grid.primToCons(rightState);
    
    ConsVar leftFlux = computePhysicalFlux(leftState);
    ConsVar rightFlux = computePhysicalFlux(rightState);

    double alpha = std::max(computeWaveSpeed(leftState), computeWaveSpeed(rightState));
    
    ConsVar flux;
    flux.rho   = 0.5 * (leftFlux.rho + rightFlux.rho)   - 0.5 * alpha * (rightCons.rho - leftCons.rho);
    flux.rhou  = 0.5 * (leftFlux.rhou + rightFlux.rhou)  - 0.5 * alpha * (rightCons.rhou - leftCons.rhou);
    flux.rhoE  = 0.5 * (leftFlux.rhoE + rightFlux.rhoE)  - 0.5 * alpha * (rightCons.rhoE - leftCons.rhoE);
    
    return flux;
}

double RusanovSolver::computeTimeStep() {
    auto primitives = grid.getPrimitives();
    double maxSpeed = 0.0;
    double dx = grid.getDx();

    for (const auto& state : primitives) {
        maxSpeed = std::max(maxSpeed, computeWaveSpeed(state));
    }
    
    return config.cfl * dx / maxSpeed;
}

void RusanovSolver::evolve() {
    double currentTime = 0.0;
    int step = 0;

    std::cout << "start t = " << config.tEnd << std::endl;

    while (currentTime < config.tEnd) {
        double dt = computeTimeStep();
        if (currentTime + dt > config.tEnd) {
            dt = config.tEnd - currentTime;
        }

        auto primitives = grid.getPrimitives();
        int nx = grid.getNx();
        std::vector<ConsVar> fluxes(nx);

        for (int i = 1; i < nx; ++i) {
            fluxes[i] = computeRusanovFlux(primitives[i-1], primitives[i]);
        }

        grid.updateSolution(fluxes, dt);
        currentTime += dt;
        step++;

        if (step % 10 == 0) {
            std::cout << "步数: " << step 
                     << " | 时间: " << std::fixed << std::setprecision(4) 
                     << currentTime << "s" << std::endl;
        }
    }

    std::cout << "推进完成！总步数：" << step << std::endl;
}