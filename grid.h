#pragma once
#include <vector>
#include <string>
#include "constants.h"
#include "config.h"

class Grid {
private:
    int nx;                     
    double dx;                  
    std::vector<double> xCoords;      
    std::vector<ConsVar> conservativeVars;    
    std::vector<ConsVar> conservativeVarsNew;
    Config config;

public:
    ConsVar primToCons(const State& prim) const;
    State consToPrim(const ConsVar& cons) const;
   

    Grid(int n, const Config& cfg);

    void initialize();
    void applyBoundaryConditions();
    
    int getNx() const;
    double getDx() const;
    std::vector<State> getPrimitives() const;
    const std::vector<double>& getXCoords() const;
    const std::vector<ConsVar>& getConservativeVars() const;
    void setConservativeVars(const std::vector<ConsVar>& vars);
    
    void updateSolution(const std::vector<ConsVar>& flux, double dt);
    void saveToFile(const std::string& filename) const;
};