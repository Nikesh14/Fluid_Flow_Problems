#include <bits/stdc++.h>
using namespace std;
typedef vector<vector<double>> arr;

void output_data(arr stream_func, arr vort_func, arr velocity_x, arr velocity_y, arr pressure, int nx, int ny, double dx, double dy){
    
    ofstream u_middle;
    u_middle.open("u_middle.txt");
    u_middle << "VARIABLES=\"y \", \" velocity_x \" \n";
    for(int i=ny-1; i>=0 ; i--){
        double pos = i*dy;
        u_middle << pos <<" "<<  velocity_x[i][nx/2] <<"\n";
    }
    u_middle.close();

    ofstream v_middle;
    v_middle.open("v_middle.txt");
    v_middle << "VARIABLES=\"x \", \" velocity_y \" \n";
    for(int i=nx-1; i>=0 ; i--){
        double pos = i*dx;
        v_middle << pos <<" "<<  velocity_y[ny/2][i] <<"\n";
    }
    v_middle.close();

    ofstream Lid_Driven_Cavity;
    Lid_Driven_Cavity.open("Lid_Driven_Cavity.txt");
    Lid_Driven_Cavity << "VARIABLES=\"X \",\"Y \",\"Stream_func \", \" vort_func \", \" velocity_x \", \" velocity_y \", \" pressure \" \n";
    for (int i = ny - 1; i >= 0; i--){
        for (int j = nx - 1; j >= 0; j--){
            double xpos, ypos;
            xpos = i * dy;
            ypos = j * dx;
            Lid_Driven_Cavity << xpos <<" "<< ypos <<" "<< stream_func[i][j] <<" "<< vort_func[i][j] <<" "<< velocity_x[i][j] <<" "<< velocity_y[i][j] <<" "<< pressure[i][j] <<"\n";
        }
    }
    Lid_Driven_Cavity.close();
}