#include "physics.hpp"

#include <algorithm>
#include <cmath>

#ifdef USE_OPENMP
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

#ifdef USE_OPENMP
#pragma omp parallel for reduction(+ : energy)
#endif
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 1; index < point_count - 1; ++index) {
            const double dx = grid.x[index + 1] - grid.x[index];
            const Vec2 derivative = (psi[curve][index + 1] - psi[curve][index]) / dx;
            energy += 0.5 * (derivative.x * derivative.x + derivative.y * derivative.y) * dx;
        }
    }

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
    const double inv_dt = 1.0 / dt;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < point_count; ++index) {
            const Vec2 temporal = (psi_new[curve][index] - psi_old[curve][index]) * inv_dt;

            const Vec2 diffusion =
                (laplacian(psi_new[curve], grid, index) + laplacian(psi_old[curve], grid, index)) *
                0.5;

            Vec2 force(0.0, 0.0);
            for (int other = 0; other < Ncurves; ++other) {
                if (other == curve) {
                    continue;
                }

                force = force + interaction(psi_new[curve][index], psi_new[other][index]) +
                        interaction(psi_old[curve][index], psi_old[other][index]);
            }

            residual[curve][index] = temporal - diffusion - force * 0.5;
        }
    }
}
