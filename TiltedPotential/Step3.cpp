#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

const double PI = 3.141592653589793;
const double theta1 = 3.67305;  // You might want to set this value
const double theta2 = 5.74858;  // You might want to set this value

int main() {
    // Load data from traj.csv
    std::ifstream file("traj.csv");
    std::string line;
    std::vector<double> theta_values;
    
    if (file.is_open()) {
        // Skip header line
        std::getline(file, line);
        
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string value;
            
            // Skip step value
            std::getline(iss, value, ',');
            // Read theta value
            std::getline(iss, value, ',');
            
            theta_values.push_back(std::stod(value));
        }
        file.close();
    } else {
        std::cerr << "Unable to open traj.csv for reading!" << std::endl;
        return 1;
    }
    
    // Compute the total difference in unwrapped angles
    double total_angle_difference = theta_values.back() - theta_values.front();
    // Calculate the number of full rotations by dividing by 2*pi
    double number_of_rotations = total_angle_difference / (2 * PI);
    std::cout << "Number of Rotations: " << number_of_rotations << std::endl;

    // Variables for processing region entries and exits
    bool in_region = false;
    double entry_time = 0.0;
    std::vector<double> time_spent_array_optimized;
    std::vector<double> exit_times;
    std::vector<double> entry_times;

    // Iterate over theta_values
    for (size_t idx = 0; idx < theta_values.size(); ++idx) {
        double theta = std::fmod(theta_values[idx], 2 * PI);
        double step = static_cast<double>(idx);  // Assuming step increments by 1

        // Check if the particle enters the region
        if (theta1 <= theta && theta <= theta1 + 5e-5 && !in_region) {
            in_region = true;
            entry_time = step;
            entry_times.push_back(entry_time);
        }
        // Check if the particle exits the region
        else if (theta > theta2 && in_region) {
            in_region = false;
            double exit_time = step;
            double time_spent = exit_time - entry_time;
            time_spent_array_optimized.push_back(time_spent);
            exit_times.push_back(exit_time);
        }
    }

    // Print the number of trajectories
    std::cout << "Number of Trajectories: " << time_spent_array_optimized.size() << std::endl;

    // Open a file to save the results
    std::ofstream outFile("results.csv");
    if (!outFile) {
        std::cerr << "Unable to open results.csv for writing!" << std::endl;
        return 1;
    }

    // Write headers
    outFile << "Entry Times,Exit Times,Time Spent\n";

    // Determine the maximum length among the vectors to iterate safely
    size_t max_length = std::max({entry_times.size(), exit_times.size(), time_spent_array_optimized.size()});

    // Write data
    for (size_t i = 0; i < max_length; ++i) {
        if (i < entry_times.size()) {
            outFile << entry_times[i];
        }
        outFile << ",";

        if (i < exit_times.size()) {
            outFile << exit_times[i];
        }
        outFile << ",";

        if (i < time_spent_array_optimized.size()) {
            outFile << time_spent_array_optimized[i];
        }
        outFile << "\n";
    }

    outFile.close();
    std::cout << "Data saved to results.csv." << std::endl;

    return 0;
}
