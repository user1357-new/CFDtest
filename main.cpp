#include <iostream>
#include <string>
#include "config.h"
#include "grid.h"
#include "rusano.h"
#include "fvs_solver.h"

int main(int argc, char* argv[]) {
    Config config;
    config.parseCommandLine(argc, argv);
    config.print();

    Grid grid(config.defaultNx, config);
    
    std::string initialFile = config.outputDir + "initial.dat";
    grid.saveToFile(initialFile);
    std::cout << "save " << initialFile << std::endl;

    // 选择求解器：1-一阶FVS(欧拉), 2-二阶FVS(欧拉), 3-一阶FVS(RK3), 4-二阶FVS(RK3), 5-Rusanov
    int solverType = 3;  // 默认使用二阶FVS + RK3
    
    if (solverType == 1) {
        std::cout << "=== 一阶 FVS + Euler ===" << std::endl;
        FVSSolver solver(grid, config, 1, 1);
        solver.evolve();
    } else if (solverType == 2) {
        std::cout << "=== 二阶 FVS + Euler ===" << std::endl;
        FVSSolver solver(grid, config, 2, 1);
        solver.evolve();
    } else if (solverType == 3) {
        std::cout << "=== 一阶 FVS + RK3 ===" << std::endl;
        FVSSolver solver(grid, config, 1, 2);
        solver.evolve();
    } else if (solverType == 4) {
        std::cout << "=== 二阶 FVS + RK3 ===" << std::endl;
        FVSSolver solver(grid, config, 2, 2);
        solver.evolve();
    } else {
        std::cout << "=== Rusanov Solver ===" << std::endl;
        RusanovSolver solver(grid, config);
        solver.evolve();
    }

    std::string resultFile = config.outputDir + "shock_tube_result_" + 
                            std::to_string(config.defaultNx) + ".dat";
    grid.saveToFile(resultFile);

    std::cout << "result saved to " << resultFile << std::endl;
    return 0;
}