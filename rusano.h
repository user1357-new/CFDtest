#pragma once
#include "grid.h"

class RusanovSolver {
private:
    Grid& grid;
    Config config;

    ConsVar computePhysicalFlux(const State& state) const;
    double computeWaveSpeed(const State& state) const;
    double computeTimeStep();

public:
    ConsVar computeRusanovFlux(const State& leftState, const State& rightState);    
    RusanovSolver(Grid& g, const Config& cfg);
    void evolve();
};