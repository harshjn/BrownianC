// In this code, I am going to look at the variation of the transition times for a tilted potential. I will feed in untilted potential and the force separately. 
// Step 0 is to generate a potential in python or externally, and find it's barrier height, maxima and minima points. 
// Step 1 of this code is to take input of a potential and then interpolate it. 
// Step 2, following that, you have to generate really long N-vector geometry for N-particles. 
// Step 3, is to use the information about maxima and minima points to generate insights back in python, by analysing the trajectories. 
//

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>

const double a = 2e-6; // size of the particle
const double R = 10e-6;  // 10 micrometers Radius of Orbit
const double eta = 1.002e-3;  // Viscosity of water at 20°C
const double wallFactor=1;
const double drag_coefficient = wallFactor*6 * M_PI * eta * a;
const double k_B = 1.38e-23;  // Boltzmann constant in J/K
const double T = 293;  // Absolute temperature (20°C in Kelvin)
const double dt = 1e-1;  // Time step
const double A = 10 * k_B * T; // 20kbT will be potential depth therefore
const double dtheta = 1e-5;
const double F0 = 2.1e-15; // Drive force
const int nums = 1e7;


int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);
    
    double theta_sim = 0.0;

    // Open traj.csv for writing
    std::ofstream trajFile("traj.csv");
    if(!trajFile) {
        std::cerr << "Unable to open traj.csv for writing!" << std::endl;
        return 1;
    }
    trajFile << "Step,Theta\n";  // header for traj.csv

    for(int i = 0; i < nums; ++i) {
        double F_theta = A * sin(theta_sim)+R*F0;
        theta_sim += (dt * F_theta / drag_coefficient + sqrt(2 * k_B * T * dt / drag_coefficient) * d(gen)) / R;
        trajFile << i << "," << theta_sim << "\n";
    }

    trajFile.close();
    std::cout << "Data saved to Pot.csv and traj.csv." << std::endl;

    return 0;
}
