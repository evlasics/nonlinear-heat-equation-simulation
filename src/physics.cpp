#include "physics.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

Vec2 interaction(const Vec2& a, const Vec2& b)
{
    const Vec2 delta = a - b;
    const double r2 = norm2(delta) + eps_reg * eps_reg;
    return delta / r2;
}

Vec2 laplacian(const std::vector<Vec2>& u, const Grid& grid, int index)
{
    const int point_count = grid.size();
    if (index == 0 || index == point_count - 1) {
        return Vec2(0.0, 0.0);
    }

    const double dx_left = grid.x[index] - grid.x[index - 1];
    const double dx_right = grid.x[index + 1] - grid.x[index];

    const Vec2 term =
        (u[index + 1] - u[index]) / dx_right - (u[index] - u[index - 1]) / dx_left;

    return term * (2.0 / (dx_left + dx_right));
}

double compute_energy(const CurveSet& psi, const Grid& grid)
{
    const int point_count = grid.size();
    double energy = 0.0;

#pragma omp parallel for reduction(+ : energy) if (point_count > 128)
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 1; index < point_count - 1; ++index) {
            const double dx = grid.x[index + 1] - grid.x[index];
            const Vec2 derivative = (psi[curve][index + 1] - psi[curve][index]) / dx;
            energy += 0.5 * (derivative.x * derivative.x + derivative.y * derivative.y) * dx;
        }
    }

#pragma omp parallel for reduction(+ : energy) if (point_count > 128)
    for (int curve_a = 0; curve_a < Ncurves; ++curve_a) {
        for (int curve_b = curve_a + 1; curve_b < Ncurves; ++curve_b) {
            for (int index = 0; index < point_count; ++index) {
                const Vec2 diff = psi[curve_a][index] - psi[curve_b][index];
                const double r2 = norm2(diff) + eps_reg * eps_reg;
                energy -= 0.5 * std::log(r2);
            }
        }
    }

    return energy;
}

void compute_residual(
    const CurveSet& psi_old,
    const CurveSet& psi_new,
    const Grid& grid,
    double dt,
    CurveSet& residual)
{
    const int point_count = grid.size();

#pragma omp parallel for if (point_count > 128)
    for (int curve = 0; curve < Ncurves; ++curve) {
        residual[curve][0] = (psi_new[curve][0] - psi_old[curve][0]) * (1.0 / dt);
        for (int index = 1; index < point_count - 1; ++index) {
            const Vec2 temporal = (psi_new[curve][index] - psi_old[curve][index]) * (1.0 / dt);
            const Vec2 diffusion =
                (laplacian(psi_new[curve], grid, index) + laplacian(psi_old[curve], grid, index)) *
                0.5;
            residual[curve][index] = temporal - diffusion;
        }
        if (point_count > 1) {
            residual[curve][point_count - 1] =
                (psi_new[curve][point_count - 1] - psi_old[curve][point_count - 1]) * (1.0 / dt);
        }
    }

#pragma omp parallel for if (point_count > 128)
    for (int index = 0; index < point_count; ++index) {
        std::vector<Vec2> force(Ncurves, Vec2(0.0, 0.0));

        for (int curve_a = 0; curve_a < Ncurves; ++curve_a) {
            for (int curve_b = curve_a + 1; curve_b < Ncurves; ++curve_b) {
                const Vec2 new_force = interaction(psi_new[curve_a][index], psi_new[curve_b][index]);
                const Vec2 old_force = interaction(psi_old[curve_a][index], psi_old[curve_b][index]);
                const Vec2 pair_force = new_force + old_force;

                force[curve_a] = force[curve_a] + pair_force;
                force[curve_b] = force[curve_b] - pair_force;
            }
        }

        for (int curve = 0; curve < Ncurves; ++curve) {
            residual[curve][index] = residual[curve][index] - force[curve] * 0.5;
        }
    }
}
