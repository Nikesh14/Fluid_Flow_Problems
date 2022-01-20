#include <bits/stdc++.h>
#include "Solver.cpp"
using namespace std;
typedef vector<vector<double>> arr;

int main(int argc, char** argv){

    //double length = 1.0, height = 1.0;
    //double dx = 0.0078125, dy = 0.0078125, dt;
    double length, height, density;
    double dx, dy, dt;
    cout<<"The length of the domain is:- "; cin>> length;
    cout<<"The height of the domain is:- "; cin>> height;
    cout<<"The density of the fluid is:- "; cin>> density;
    cout<<"Spatial discretization in x:- "; cin>> dx;
    cout<<"Spatial discretization in y:- "; cin>> dy;
    cout<<"Spatial discretization in time:- "; cin>> dt;
    double convergence; cout<<"The minimum residual error:- "; cin>> convergence ; 
    double reynolds_number; cout<<"Reynolds number is:- "; cin>> reynolds_number; 
    double relaxtion; cout<<"relaxation factor = 1, Gauss Seidel Method \n" 
                          <<"relaxation factor < 1, Under-relaxation Method \n" 
                          <<"relaxation factor > 1, Over-relaxation Method \n"
                          <<"Enter the relaxation factor:- "; cin>> relaxtion; 
    double nx = (length / dx) + 1, ny = (height / dy) + 1;
    arr velocity_x, velocity_y;
    arr stream_func, vort_func;
    arr pressure;

    // initializing Velocity fields---------------------
    for (int i = 0; i < ny; ++i){
        vector<double> temp_x, temp_y;
        for (int j = 0; j < nx; ++j){
            if (i == 0)
                temp_x.push_back(1);
            else
                temp_x.push_back(0);
            temp_y.push_back(0);
        }
        velocity_x.push_back(temp_x);
        velocity_y.push_back(temp_y);
        vort_func.push_back(temp_y);
        stream_func.push_back(temp_y);
        pressure.push_back(temp_y);
    }
    Solve_Prb(velocity_x, velocity_y, vort_func, stream_func, pressure, nx, ny, dx, dy, dt, density, convergence, reynolds_number, relaxtion);
    return 0;
}