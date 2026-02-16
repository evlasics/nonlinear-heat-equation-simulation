#include <algorithm>
#include <iostream>

#include "config.hpp"
#include "output.hpp"
#include "simulation_types.hpp"
#include "solver.hpp"

namespace {

constexpr int kOutputEvery = 100;
constexpr int kLogEvery = 1000;
constexpr int kRefineEvery = 10;
constexpr double kRefineThreshold = 0.5;

} // namespace

int main()
{
    Grid grid;
    grid.initialize(Nx, L);

    CurveSet psi(Ncurves, std::vector<Vec2>(Nx));
    initialize_curves(psi, grid);

    std::string out_dir = make_run_directory();
    PVDWriter pvd(out_dir);

    double dt = 1e-3;
    double time = 0.0;
    int frame = 0;

    while (time < Tfinal) {
        const TimeStepResult step = attempt_timestep(psi, grid, dt);
        if (!step.accepted) {
            if (dt < dt_min) {
                std::cout << "Time step underflow.\n";
                break;
            }
            continue;
        }

        time += dt;

        if (frame % kOutputEvery == 0) {
            write_vtp(psi, grid, frame, out_dir);
            pvd.append(frame, time);

            if (frame % kLogEvery == 0) {
                const double energy = compute_energy(psi, grid);
                std::cout << "t=" << time << "  dt=" << dt << "  Nx=" << grid.size()
                          << "  Energy=" << energy << "\n";
            }
        }

        if (step.error < tol_time * 0.1) {
            dt = std::min(dt * 1.2, dt_max);
        }

        ++frame;
        if (frame % kRefineEvery == 0) {
            refine_mesh(grid, psi, kRefineThreshold);
        }
    }

    pvd.finalize();
    archive_latest_run();
    std::cout << "Finished.\n";

    return 0;
}
