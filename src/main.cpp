#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <filesystem>
#include "config.hpp"

// =====================================================
// PARAMETERS
// =====================================================

double dt = 1e-3;
std::string out_dir;


// =====================================================
// INTERACTION
// =====================================================

Vec2 interaction(const Vec2& a, const Vec2& b){
    Vec2 d = a-b;
    double r2 = norm2(d) + eps_reg*eps_reg;
    return d / r2;
}

// =====================================================
// ENERGY
// =====================================================

double compute_energy(
    const vector<vector<Vec2>>& psi,
    const Grid& grid)
{
    int Nx = grid.size();
    double E = 0.0;

    for(int j=0;j<Ncurves;j++){
        for(int i=1;i<Nx-1;i++){
            double dx = grid.x[i+1]-grid.x[i];
            Vec2 d = (psi[j][i+1]-psi[j][i]) / dx;
            E += 0.5*(d.x*d.x+d.y*d.y)*dx;
        }
    }

    for(int j=0;j<Ncurves;j++)
    for(int k=j+1;k<Ncurves;k++)
    for(int i=0;i<Nx;i++){
        Vec2 diff = psi[j][i]-psi[k][i];
        double r2 = norm2(diff)+eps_reg*eps_reg;
        E -= 0.5*log(r2);
    }

  return E;
}

// =====================================================
// DISCRETE LAPLACIAN (nonuniform grid)
// =====================================================

Vec2 laplacian(
    const vector<Vec2>& u,
    const Grid& grid,
    int i)
{
    int Nx = grid.size();
    if(i==0 || i==Nx-1) return Vec2(0,0);

    double dxL = grid.x[i]-grid.x[i-1];
    double dxR = grid.x[i+1]-grid.x[i];

    Vec2 term =
        (u[i+1]-u[i])/dxR
      - (u[i]-u[i-1])/dxL;

    return term * (2.0/(dxL+dxR));
}

// =====================================================
// RESIDUAL for Crank-Nicolson
// =====================================================

void compute_residual(
    const vector<vector<Vec2>>& psi_old,
    const vector<vector<Vec2>>& psi_new,
    const Grid& grid,
    vector<vector<Vec2>>& R)
{
    int Nx = grid.size();

    for(int j=0;j<Ncurves;j++)
    for(int i=0;i<Nx;i++){

        Vec2 diff = (psi_new[j][i]-psi_old[j][i])*(1.0/dt);

        Vec2 lap = (laplacian(psi_new[j],grid,i)
                   +laplacian(psi_old[j],grid,i))*0.5;

        Vec2 force(0,0);
        for(int k=0;k<Ncurves;k++){
            if(k==j) continue;
            force = force + interaction(psi_new[j][i],psi_new[k][i])
                            + interaction(psi_old[j][i],psi_old[k][i]);
        }
        force = force*0.5;

        R[j][i] = diff - lap - force;
    }
}

// =====================================================
// SIMPLE NEWTON SOLVER
// (Jacobian approximated diagonally for robustness)
// =====================================================

bool newton_step(
    vector<vector<Vec2>>& psi,
    const vector<vector<Vec2>>& psi_old,
    const Grid& grid)
{
    int Nx = grid.size();

    vector<vector<Vec2>> R(Ncurves, vector<Vec2>(Nx));

    for(int iter=0; iter<max_newton; iter++){

        compute_residual(psi_old,psi,grid,R);

        double maxR = 0.0;
        for(int j=0;j<Ncurves;j++)
        for(int i=0;i<Nx;i++){
            maxR = max(maxR, sqrt(norm2(R[j][i])));
        }

        if(maxR < tol_newton) return true;

        // Diagonal update (robust but not fast)
        for(int j=0;j<Ncurves;j++)
        for(int i=0;i<Nx;i++){
            psi[j][i] = psi[j][i] - R[j][i]*dt;
        }
    }

    return false;
}

// =====================================================
// ADAPTIVE MESH REFINEMENT
// =====================================================

void refine_mesh(
    Grid& grid,
    vector<vector<Vec2>>& psi)
{
    int Nx = grid.size();
    vector<double> new_x;
    vector<vector<Vec2>> new_psi(Ncurves);

    for(int j=0;j<Ncurves;j++)
        new_psi[j].clear();

    for(int i=0;i<Nx-1;i++){

        new_x.push_back(grid.x[i]);
        for(int j=0;j<Ncurves;j++)
            new_psi[j].push_back(psi[j][i]);

        double indicator = 0.0;

        for(int j=0;j<Ncurves;j++){
            Vec2 d = psi[j][i+1]-psi[j][i];
            indicator += norm2(d);
        }

        if(indicator > 0.5){  // refinement threshold

            double xm = 0.5*(grid.x[i]+grid.x[i+1]);
            new_x.push_back(xm);

            for(int j=0;j<Ncurves;j++){
                Vec2 mid = (psi[j][i]+psi[j][i+1])*0.5;
                new_psi[j].push_back(mid);
            }
        }
    }

    new_x.push_back(grid.x.back());
    for(int j=0;j<Ncurves;j++)
        new_psi[j].push_back(psi[j].back());

    grid.x = new_x;
    psi = new_psi;
}

// =====================================================
// CREATE OUTPUT DIRECTORY
// =====================================================

std::string make_run_directory() {
    using namespace std::chrono;

    auto now = system_clock::now();
    auto t = system_clock::to_time_t(now);

    std::stringstream ss;
    ss << "output/run_" << std::put_time(std::localtime(&t), "%Y_%m_%d_%H_%M_%S") << "/";

    std::filesystem::create_directories(ss.str());

    return ss.str();
}

// =====================================================
// WRITE VTP FILE
// =====================================================

void write_vtp(
    const vector<vector<Vec2>>& psi,
    const Grid& grid,
    int frame)
{
    int Nx = grid.size();
    int total_points = Ncurves * Nx;

    char filename[256];
		std::string path = out_dir + "frame_%04d.vtp";
    sprintf(filename, path.c_str(), frame);

    std::ofstream file(filename);

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <PolyData>\n";
    file << "    <Piece NumberOfPoints=\""
         << total_points
         << "\" NumberOfLines=\""
         << Ncurves
         << "\">\n";

    // ---- Points ----
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";

    for(int j=0; j<Ncurves; j++)
    for(int i=0; i<Nx; i++)
    {
        double x = grid.x[i];
        file << x << " "
             << psi[j][i].x << " "
             << psi[j][i].y << "\n";
    }

    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // ---- Lines connectivity ----
    file << "      <Lines>\n";

    // connectivity
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for(int j=0; j<Ncurves; j++)
    {
        for(int i=0; i<Nx; i++)
            file << j*Nx + i << " ";
        file << "\n";
    }
    file << "        </DataArray>\n";

    // offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for(int j=0; j<Ncurves; j++)
        file << (j+1)*Nx << " ";
    file << "\n";
    file << "        </DataArray>\n";

    file << "      </Lines>\n";

    file << "    </Piece>\n";
    file << "  </PolyData>\n";
    file << "</VTKFile>\n";

    file.close();
}

class PVDWriter
{
    std::ofstream file;
    bool initialized = false;

public:
		void initialize()
		{
				file.open(out_dir + "solution.pvd");
				file << "<?xml version=\"1.0\"?>\n";
				file << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
				file << "  <Collection>\n";
				initialized = true;
		}

		void append(int frame, double time)
		{
				if(!initialized) initialize();

				file << "    <DataSet timestep=\""
						 << std::setprecision(16)
						 << time
						 << "\" group=\"\" part=\"0\" file=\"frame_"
						 << std::setw(4) << std::setfill('0') << frame
						 << ".vtp\"/>\n";
		}



    void finalize()
    {
        file << "  </Collection>\n";
        file << "</VTKFile>\n";
        file.close();
    }
};



// =====================================================
// MAIN
// =====================================================

int main(){
		PVDWriter pvd;
		int frame = 0;

    Grid grid;
    grid.initialize(Nx, L);

    vector<vector<Vec2>> psi(Ncurves,
                             vector<Vec2>(Nx));

		initialize_curves(psi, grid);

		out_dir = make_run_directory();
		

    double t = 0.0;

    while(t < Tfinal){
        auto psi_old = psi;
        auto psi_trial = psi;

        bool ok = newton_step(psi_trial,
                              psi_old,
                              grid);

        if(!ok){
            dt *= 0.5;
            if(dt < dt_min){
                cout<<"Time step underflow.\n";
                break;
            }
            continue;
        }

        // adaptive timestep using difference
        double err = 0.0;
        for(int j=0;j<Ncurves;j++)
        for(int i=0;i<grid.size();i++)
            err = max(err,
                      sqrt(norm2(psi_trial[j][i]
                           -psi_old[j][i])));

        if(err > tol_time){
            dt *= 0.5;
            continue;
        }

        psi = psi_trial;
        t += dt;

				if(frame % 100 == 0)
				{
					write_vtp(psi, grid, frame);
					pvd.append(frame, t);

					if (frame % 1000 == 0)
					{
						double E = compute_energy(psi,grid);
		
						cout<<"t="<<t
								<<"  dt="<<dt
								<<"  Nx="<<grid.size()
								<<"  Energy="<<E
								<<endl;
					}
				}
				frame++;


        if(err < tol_time*0.1)
            dt = min(dt*1.2, dt_max);

				if(frame % 10 == 0)
			    refine_mesh(grid, psi);

        refine_mesh(grid,psi);

    }

		pvd.finalize();

    cout<<"Finished.\n";
}

