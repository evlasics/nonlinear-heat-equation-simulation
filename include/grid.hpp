#pragma once
#include <vector>
using std::vector;

using namespace std;

struct Grid {
	std::vector<double> x;

    void initialize(int Nx, double L){
        x.resize(Nx);
        double dx = 2.0*L/(Nx-1);
        for(int i=0;i<Nx;i++)
            x[i] = -L + i*dx;
    }

    int size() const { return x.size(); }
};
