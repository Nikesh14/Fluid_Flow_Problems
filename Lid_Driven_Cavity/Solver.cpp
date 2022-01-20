#include <bits/stdc++.h>
#include "Processing_data.cpp"
using namespace std;
typedef vector<vector<double>> arr;

void Solve_Prb(arr velocity_x, arr velocity_y, arr vort_func, arr stream_func, arr pressure, int nx, int ny, double dx, double dy, double dt, double density, double convergence, double re, double beta){

    double MAX_error = INT_MAX;
    arr stream_func_old, vort_func_old;
    stream_func_old.assign(stream_func.begin(), stream_func.end());
    vort_func_old.assign(vort_func.begin(), vort_func.end());
    int iter = 0;
    while (MAX_error>convergence){

        // Calculation of Stream Function //
        for (int i = 1; i < ny - 1; ++i){
            for (int j = 1; j < nx - 1; ++j){

                stream_func[i][j] = (0.5*beta*(dx*dx*dy*dy)/(dx*dx + dy*dy))*(
                (stream_func[i][j+1] + stream_func[i][j-1])/(dx*dx)+
                (stream_func[i+1][j] + stream_func[i-1][j])/(dy*dy) + vort_func[i][j])
                +(1-beta)*stream_func[i][j];

            }
        }

        // Calculation of Vorticity Function //
        for (int i = 1; i < ny - 1; ++i){
            for (int j = 1; j < nx - 1; ++j){
                vort_func[i][j] = vort_func[i][j]
                -dt*(((stream_func[i+1][j] - stream_func[i-1][j])/(2*dy))*((vort_func[i][j+1] - vort_func[i][j-1])/(2*dx)))
                +dt*(((stream_func[i][j+1] - stream_func[i][j-1])/(2*dx))*((vort_func[i+1][j] - vort_func[i-1][j])/(2*dy)))
                +dt*((1/re)*((vort_func[i][j+1] - 2*vort_func[i][j] + vort_func[i][j-1])/(dx*dx) + (vort_func[i+1][j] - 2*vort_func[i][j] + vort_func[i-1][j])/(dy*dy)));
            }
        }

        // Boundary Conditions //
        for (int i = 1; i < nx - 1; ++i){
            vort_func[ny-1][i] = (-2/(dy*dy))*(stream_func[ny-2][i] - stream_func[ny-1][i]);                //Bottom wall
            vort_func[0][i] = (-2/(dy*dy))*(stream_func[1][i] - stream_func[0][i]) - (2*velocity_x[0][i])/dy; //top wall
        }
        for (int j = 0; j < ny; ++j){
            vort_func[j][0] = (-2/(dx*dx))*(stream_func[j][1] - stream_func[j][0]);                //left wall
            vort_func[j][nx-1] = (-2/(dx*dx))*(stream_func[j][nx-2] - stream_func[j][nx-1]); //right wall
        }

        // Calculation of Velocities //
        for (int i = 1; i < ny-1; ++i){
            for (int j = 1; j < nx-1; ++j){
                velocity_x[i][j] = -(stream_func[i + 1][j] - stream_func[i - 1][j]) / (2 * dy);
                velocity_y[i][j] = -(stream_func[i][j + 1] - stream_func[i][j - 1]) / (2 * dx);
            }
        }

        // Calculation of Pressure //
        for (int i = 1; i < ny-1; ++i){
            for (int j = 1; j < nx-1; ++j){
                pressure[i][j] = 0.5*(dx*dx*dy*dy)/(dx*dx + dy*dy)*(
                (pressure[i][j+1] + pressure[i][j-1])/(dx*dx) + (pressure[i+1][j] + pressure[i-1][j])/(dy*dy)- 2*density*(
                ((velocity_x[i][j+1] - velocity_x[i][j-1])/(2*dx))*((velocity_y[i+1][j] - velocity_y[i-1][j])/(2*dy))-
                ((velocity_y[i][j+1] - velocity_y[i][j-1])/(2*dx))*((velocity_x[i+1][j] - velocity_x[i-1][j])/(2*dy))));
            }
        }

        // Pressure Boundary Conditions //
        for(int i=0; i<nx; ++i)
            pressure[0][i] = pressure[1][i]; // Bottom Wall
        for(int i=0; i<nx; ++i)
            pressure[ny-1][i] = pressure[ny-1][i]; // Top Wall
        
        for(int i=0; i<ny; ++i)
            pressure[i][0] = pressure[i][1]; // Left Wall
        for(int i=0; i<ny; ++i)
            pressure[i][nx-1] = pressure[i][nx-1]; // Right Wall

    
        // Error Calculations // 
        vector<double> error;
        double err = 0;
        for (int i = 0; i < ny; ++i){
            vector<double> temp;
            for (int j = 0; j < nx; ++j){
                double var = abs(abs(stream_func[i][j]) - abs(stream_func_old[i][j]));
                temp.push_back(var);
            }
            err = *max_element(temp.begin(), temp.end());
            error.push_back(err);
        }
        MAX_error = *max_element(error.begin(), error.end());
        if (MAX_error == 0)
            MAX_error = 0.1;
        cout << "Error at " << iter << " is " << MAX_error << "\n";
        ++iter;
        stream_func_old.assign(stream_func.begin(), stream_func.end());
        vort_func_old.assign(vort_func.begin(), vort_func.end());
    }
    output_data(stream_func, vort_func, velocity_x, velocity_y, pressure, nx, ny, dx, dy);
}