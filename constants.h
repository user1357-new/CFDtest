#pragma once
#include <cmath>

// 状态向量（原始变量：rho, u, p）
struct State {
    double rho, u, p;
    State() : rho(0), u(0), p(0) {}
    State(double r, double uu, double pr) : rho(r), u(uu), p(pr) {}
    
    State operator+(const State& s) const { return State(rho+s.rho, u+s.u, p+s.p); }
    State operator-(const State& s) const { return State(rho-s.rho, u-s.u, p-s.p); }
    State operator*(double c) const { return State(rho*c, u*c, p*c); }
    State operator/(double c) const { return State(rho/c, u/c, p/c); }
};

// 保守变量向量（rho, rho*u, rho*E）
struct ConsVar {
    double rho, rhou, rhoE;
    ConsVar() : rho(0), rhou(0), rhoE(0) {}
    ConsVar(double r, double ru, double re) : rho(r), rhou(ru), rhoE(re) {}

    ConsVar operator+(const ConsVar& s) const { return ConsVar(rho+s.rho, rhou+s.rhou, rhoE+s.rhoE); }
    ConsVar operator-(const ConsVar& s) const { return ConsVar(rho-s.rho, rhou-s.rhou, rhoE-s.rhoE); }
    ConsVar operator*(double c) const { return ConsVar(rho*c, rhou*c, rhoE*c); }
    ConsVar operator/(double c) const { return ConsVar(rho/c, rhou/c, rhoE/c); }
};