// config.cpp
#include "config.hpp"

#include <cmath>

double norm2(const Vec2& v){ return v.x*v.x+v.y*v.y; }

int Ncurves = 2;
int Nx = 400;
double L = 20.0;
double Tfinal = 2.0;

double dt_min = 1e-6;
double dt_max = 1e-1;

double eps_reg = 1e-8;
double tol_newton = 1e-8;
int max_newton = 20;
double tol_time = 1e-5;

void initialize_curves(
    std::vector<std::vector<Vec2>>& psi,
    const Grid& grid)
{
    for(int i=0; i<Nx; i++)
    {
        double x = grid.x[i];

        psi[0][i] = Vec2(0.1*std::cos(x), x + 1);
        psi[1][i] = Vec2(x, 0.1*std::sin(x));
    }
}
