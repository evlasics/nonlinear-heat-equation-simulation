// config.h
#pragma once
#include "vector.hpp"
#include "grid.hpp"

extern int Nx;
extern int Ncurves;

extern double output_dt;

extern double L;
extern double Tfinal;

extern double dt_min;
extern double dt_max;

extern double eps_reg;
extern double tol_newton;
extern int max_newton;
extern double tol_time;

void initialize_curves(
    std::vector<std::vector<Vec2>>& psi,
    const Grid& grid);

double norm2(const Vec2& v);
