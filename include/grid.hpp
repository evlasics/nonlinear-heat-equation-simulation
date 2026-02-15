#pragma once

#include <vector>

struct Grid {
    std::vector<double> x;

    void initialize(int Nx, double L)
    {
        x.resize(Nx);
        const double dx = 2.0 * L / (Nx - 1);

        for (int i = 0; i < Nx; i++) {
            x[i] = -L + i * dx;
        }
    }

    int size() const { return static_cast<int>(x.size()); }
};
