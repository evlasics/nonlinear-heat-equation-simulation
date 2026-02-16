#include "output.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <utility>

namespace {

std::string make_archive_directory()
{
    using namespace std::chrono;

    const auto now = system_clock::now();
    const auto time_now = system_clock::to_time_t(now);
    const auto millis = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    std::ostringstream name;
    name << "output/archive/run_" << std::put_time(std::localtime(&time_now), "%Y_%m_%d_%H_%M_%S")
         << '_' << std::setw(3) << std::setfill('0') << millis.count();

    std::filesystem::create_directories(name.str());
    return name.str() + "/";
}

} // namespace

void clear_run_outputs(const std::string& dir)
{
    if (!std::filesystem::exists(dir)) {
        return;
    }

    for (const auto& entry : std::filesystem::directory_iterator(dir)) {
        if (!entry.is_regular_file()) {
            continue;
        }

        const std::string name = entry.path().filename().string();
        const bool is_frame =
            name.rfind("frame_", 0) == 0 && entry.path().extension() == ".vtp";

        if (is_frame || name == "solution.pvd") {
            std::filesystem::remove(entry.path());
        }
    }
}

std::string make_run_directory()
{
    const std::string dir = "output/latest/";
    std::filesystem::create_directories(dir);
    clear_run_outputs(dir);
    return dir;
}

void archive_latest_run()
{
    const std::filesystem::path latest_dir("output/latest");
    if (!std::filesystem::exists(latest_dir)) {
        return;
    }

    const std::filesystem::path archive_dir = make_archive_directory();
    std::filesystem::copy(
        latest_dir,
        archive_dir,
        std::filesystem::copy_options::recursive | std::filesystem::copy_options::overwrite_existing);
}

void write_vtp(const CurveSet& psi, const Grid& grid, int frame, const std::string& out_dir)
{
    const int point_count = grid.size();
    const int total_points = Ncurves * point_count;

    CurveSet second_derivative(Ncurves, std::vector<Vec2>(point_count, Vec2(0.0, 0.0)));
    std::vector<std::vector<double>> distance_to_other_curve(
        Ncurves,
        std::vector<double>(point_count, 0.0));

    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < point_count; ++index) {
            second_derivative[curve][index] = laplacian(psi[curve], grid, index);

            double min_dist = std::numeric_limits<double>::max();
            for (int other = 0; other < Ncurves; ++other) {
                if (other == curve) {
                    continue;
                }

                const Vec2 delta = psi[curve][index] - psi[other][index];
                min_dist = std::min(min_dist, std::sqrt(norm2(delta)));
            }

            if (min_dist < std::numeric_limits<double>::max()) {
                distance_to_other_curve[curve][index] = min_dist;
            }
        }
    }

    char file_name[256];
    const std::string path_pattern = out_dir + "frame_%04d.vtp";
    std::snprintf(file_name, sizeof(file_name), path_pattern.c_str(), frame);

    std::ofstream file(file_name);

    // XML header.
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <PolyData>\n";
    file << "    <Piece NumberOfPoints=\"" << total_points << "\" NumberOfLines=\"" << Ncurves
         << "\">\n";

    // 3D points: x-grid coordinate and 2D curve value in y/z slots.
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < point_count; ++index) {
            file << grid.x[index] << ' ' << psi[curve][index].x << ' ' << psi[curve][index].y << "\n";
        }
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // Polyline connectivity for each curve.
    file << "      <Lines>\n";
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < point_count; ++index) {
            file << curve * point_count + index << ' ';
        }
        file << "\n";
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int curve = 0; curve < Ncurves; ++curve) {
        file << (curve + 1) * point_count << ' ';
    }
    file << "\n";
    file << "        </DataArray>\n";
    file << "      </Lines>\n";

    // Derived fields for visualization and diagnostics.
    file << "      <PointData Scalars=\"distance_to_other_curve second_derivative_magnitude\">\n";

    file << "        <DataArray type=\"Float64\" Name=\"second_derivative\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < point_count; ++index) {
            const Vec2& d2 = second_derivative[curve][index];
            file << d2.x << ' ' << d2.y << " 0\n";
        }
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"second_derivative_magnitude\" format=\"ascii\">\n";
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < point_count; ++index) {
            file << std::sqrt(norm2(second_derivative[curve][index])) << "\n";
        }
    }
    file << "        </DataArray>\n";

    file << "        <DataArray type=\"Float64\" Name=\"distance_to_other_curve\" format=\"ascii\">\n";
    for (int curve = 0; curve < Ncurves; ++curve) {
        for (int index = 0; index < point_count; ++index) {
            file << distance_to_other_curve[curve][index] << "\n";
        }
    }
    file << "        </DataArray>\n";

    file << "      </PointData>\n";
    file << "    </Piece>\n";
    file << "  </PolyData>\n";
    file << "</VTKFile>\n";
}

PVDWriter::PVDWriter(std::string out_dir)
    : out_dir_(std::move(out_dir)), initialized_(false)
{
}

void PVDWriter::initialize()
{
    file_.open(out_dir_ + "solution.pvd");
    file_ << "<?xml version=\"1.0\"?>\n";
    file_ << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file_ << "  <Collection>\n";
    initialized_ = true;
}

void PVDWriter::append(int frame, double time)
{
    if (!initialized_) {
        initialize();
    }

    file_ << "    <DataSet timestep=\"" << std::setprecision(16) << time
           << "\" group=\"\" part=\"0\" file=\"frame_" << std::setw(4) << std::setfill('0')
           << frame << ".vtp\"/>\n";
}

void PVDWriter::finalize()
{
    if (!initialized_) {
        return;
    }

    file_ << "  </Collection>\n";
    file_ << "</VTKFile>\n";
    file_.close();
}
