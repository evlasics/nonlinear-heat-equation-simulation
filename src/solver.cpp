#include "solver.hpp"

#include <algorithm>
#include <cmath>

namespace {

double max_residual_norm(const CurveSet& residual)
{
    double max_norm = 0.0;

    for (const auto& curve : residual) {
        for (const auto& value : curve) {
            max_norm = std::max(max_norm, std::sqrt(norm2(value)));
        }
    }

    return max_norm;
}

double compute_step_error(const CurveSet& next, const CurveSet& previous)
{
    double error = 0.0;

    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < static_cast<int>(next[curve].size()); ++index) {
            error = std::max(error, std::sqrt(norm2(next[curve][index] - previous[curve][index])));
        }
    }

    return error;
}

} // namespace

bool newton_step(CurveSet& psi_trial, const CurveSet& psi_old, const Grid& grid, double dt)
{
    const int point_count = grid.size();
    CurveSet residual(Ncurves, std::vector<Vec2>(point_count));

    for (int iteration = 0; iteration < max_newton; ++iteration) {
        compute_residual(psi_old, psi_trial, grid, dt, residual);

        if (max_residual_norm(residual) < tol_newton) {
            return true;
        }

        // Robust diagonal-only correction.
        for (int curve = 0; curve < Ncurves; ++curve) {
            for (int index = 0; index < point_count; ++index) {
                psi_trial[curve][index] = psi_trial[curve][index] - residual[curve][index] * dt;
            }
        }
    }

    return false;
}

void refine_mesh(Grid& grid, CurveSet& psi, double threshold)
{
    const int point_count = grid.size();
    std::vector<double> refined_x;
    CurveSet refined_psi(Ncurves);

    for (int index = 0; index < point_count - 1; ++index) {
        refined_x.push_back(grid.x[index]);
        for (int curve = 0; curve < Ncurves; ++curve) {
            refined_psi[curve].push_back(psi[curve][index]);
        }

        double indicator = 0.0;
        for (int curve = 0; curve < Ncurves; ++curve) {
            const Vec2 diff = psi[curve][index + 1] - psi[curve][index];
            indicator += norm2(diff);
        }

        if (indicator > threshold) {
            const double mid = 0.5 * (grid.x[index] + grid.x[index + 1]);
            refined_x.push_back(mid);

            for (int curve = 0; curve < Ncurves; ++curve) {
                const Vec2 midpoint_value = (psi[curve][index] + psi[curve][index + 1]) * 0.5;
                refined_psi[curve].push_back(midpoint_value);
            }
        }
    }

    refined_x.push_back(grid.x.back());
    for (int curve = 0; curve < Ncurves; ++curve) {
        refined_psi[curve].push_back(psi[curve].back());
    }

    grid.x = std::move(refined_x);
    psi = std::move(refined_psi);
}

TimeStepResult attempt_timestep(CurveSet& psi, Grid& grid, double& dt)
{
    const CurveSet psi_old = psi;
    CurveSet psi_trial = psi;

    if (!newton_step(psi_trial, psi_old, grid, dt)) {
        dt *= 0.5;
        return {.accepted = false, .converged = false, .error = 0.0};
    }

    const double error = compute_step_error(psi_trial, psi_old);
    if (error > tol_time) {
        dt *= 0.5;
        return {.accepted = false, .converged = true, .error = error};
    }

    psi = std::move(psi_trial);
    return {.accepted = true, .converged = true, .error = error};
}
