#include "config.h"
#include <iostream>
#include <cmath>

Config::Config() 
    : gamma(1.4),
      xMin(-0.5),
      xMax(0.5),
      tEnd(0.3),
      cfl(0.1),
      mu(0.01),
      mu4(0.001),
      defaultNx(200),
      outputDir("./") {
}

void Config::parseCommandLine(int argc, char* argv[]) {
    if (argc > 1) {
        defaultNx = std::stoi(argv[1]);
    }
}

void Config::print() const {
    std::cout << "Config" << std::endl;
    std::cout << "zone[" << xMin << ", " << xMax << "]" << std::endl;
    std::cout << "mesh " << defaultNx << std::endl;
    std::cout << "ending time" << tEnd << std::endl;
    std::cout << "CFL " << cfl << std::endl;
    std::cout << "Gamma " << gamma << std::endl;
    std::cout << "================" << std::endl;
}