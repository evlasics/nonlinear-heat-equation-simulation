#pragma once

#include <fstream>
#include <string>

#include "physics.hpp"
#include "simulation_types.hpp"

// Remove transient frame outputs from a run directory.
void clear_run_outputs(const std::string& dir);

// Prepare and clean output/latest for a new run.
std::string make_run_directory();

// Copy output/latest to output/archive with a timestamped folder.
void archive_latest_run();

// Write a frame to VTK PolyData format.
void write_vtp(const CurveSet& psi, const Grid& grid, int frame, const std::string& out_dir);

class PVDWriter {
public:
    explicit PVDWriter(std::string out_dir);
    void append(int frame, double time);
    void finalize();

private:
    void initialize();

    std::string out_dir_;
    std::ofstream file_;
    bool initialized_;
};
