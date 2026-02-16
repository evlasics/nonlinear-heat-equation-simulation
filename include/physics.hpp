#pragma once

#include "config.hpp"
#include "simulation_types.hpp"

// Pairwise regularized interaction between curve points.
Vec2 interaction(const Vec2& a, const Vec2& b);

// Discrete Laplacian on a non-uniform 1D grid.
Vec2 laplacian(const std::vector<Vec2>& u, const Grid& grid, int index);

// Total energy of the current curve configuration.
double compute_energy(const CurveSet& psi, const Grid& grid);

// Crank-Nicolson residual used by the nonlinear solver.
void compute_residual(
    const CurveSet& psi_old,
    const CurveSet& psi_new,
    const Grid& grid,
    double dt,
    CurveSet& residual);
