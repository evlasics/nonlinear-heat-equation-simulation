#pragma once

#include "physics.hpp"
#include "simulation_types.hpp"

struct TimeStepResult {
    bool accepted;
    bool converged;
    double error;
};

// Perform nonlinear solve for one implicit time step with diagonal Newton updates.
bool newton_step(CurveSet& psi_trial, const CurveSet& psi_old, const Grid& grid, double dt);

// Refine intervals whose aggregate curve jump exceeds the threshold.
void refine_mesh(Grid& grid, CurveSet& psi, double threshold);

// Attempt to advance one time step and return status info for adaptation logic.
TimeStepResult attempt_timestep(CurveSet& psi, Grid& grid, double& dt);
