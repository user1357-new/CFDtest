// 新文件: fvs_solver.h
#pragma once
#include "grid.h"
#include <vector>

class FVSSolver {
private:
    Grid& grid;
    Config config;
    int order;          // 空间精度: 1-一阶, 2-二阶
    int timeScheme;     // 时间格式: 1-欧拉, 2-RK3
    
    ConsVar computeStegerWarmingFlux(const State& state) const;
    ConsVar computeRusanovFlux(const State& left, const State& right) const;
    double computeWaveSpeed(const State& state) const;
    
    // 二阶精度需要的限制器
    double minmod(double a, double b) const;
    State minmod(const State& a, const State& b) const;
    State reconstructState(int i, const std::vector<State>& primitives, bool leftSide) const;
    
    // 时间推进方法
    void eulerStep(std::vector<ConsVar>& U_new, const std::vector<ConsVar>& U_old, double dt);
    void rk3Step(std::vector<ConsVar>& U_final, const std::vector<ConsVar>& U_old, double dt);

public:
    FVSSolver(Grid& g, const Config& cfg, int spatialOrder = 2, int timeMethod = 2);
    void evolve();
    double computeTimeStep();
};