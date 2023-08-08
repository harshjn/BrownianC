#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
//position of the target, given (50,50) is starting location of walker
// structure to represent the state of the walker


struct Walker {
    double x, y;
};
struct Target {
    double x,y;
};

int N=1000;  //Lattice Dimension
const double a = 1.0;  //Target size dimension
Target t = {N/2+13,N/2}; //Target location
// number of simulations
const int numSimulations = 10;
const int maxSteps = 1000000;
const float beta = 1;
 
Eigen::SparseMatrix<double> G(N, N);  // initialize a 100x100 sparse matrix

double distanceToTarget(Walker w, Target t) {
    double dx = w.x - t.x;
    double dy = w.y - t.y;
    return std::sqrt(dx * dx + dy * dy);
}


//Prepare random number generation
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> uniformDist(0, 1);



Walker stepDir(Walker w, Eigen::Matrix3d localWeights) {
    // replace w by e^[beta*w] in the weight matrix.
    Eigen::Matrix3d funcWeights = localWeights.unaryExpr([](double x) { return std::exp(-1*beta * x); });
    // a,b,c,d,e,f,g,h,i 9 states and 9 weights
    
    std::vector<int> directions = {-1, 0, 1};
    std::vector<double> weightsList;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            weightsList.push_back(funcWeights(i,j));

// Generate a distribution based on the weights
    std::discrete_distribution<> dist(weightsList.begin(), weightsList.end());

    // Generate a random index
    int randomIndex = dist(gen);

    // Calculate new coordinates
    int dirX = randomIndex / 3 - 1, dirY = randomIndex % 3 - 1;
 //   std::cout << "dirX: " << dirX << ", dirY: " << dirY << '\n';
   
    double stepX =uniformDist(gen);
    double stepY =uniformDist(gen);
    
    double walkerXOut = w.x + stepX*dirX;
    double walkerYOut = w.y + stepY*dirY;

    Walker walkerOut ={walkerXOut,walkerYOut};

    return walkerOut;
}

Walker jumpCalculate(Walker wIn) {
    int startX = (int) wIn.x;
    int startY = (int) wIn.y;

    Eigen::Matrix3d localWeights;
 
    for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {

        int XX = startX-1+i; 
        int YY= startY-1+j; 
        if (XX<0){
            XX+=N;
        }
        if (XX>(N-1)){
            XX-=N;
        }if (YY<0){
            YY+=N;
        }if (YY>(N-1)){
            YY-=N;
        }
        localWeights(i,j) = (float)G.coeff(XX,YY);
    }
    }
    Walker walkerOut = stepDir(wIn,localWeights);
    // std::cout << "Matrix:\n" << localWeights << std::endl;
    
    // std::cout<< startX<<std::endl;
    // Now, we pass the localWeights and position to 

    // Update G matrix
    return walkerOut;
}

// Function to write sparse matrix to text file
void writeSparseTxt(const Eigen::SparseMatrix<double>& G, int simulationNumber) {
    std::string filename = "Gsim" + std::to_string(simulationNumber) + ".txt";
    std::ofstream out(filename);
    for (int k=0; k < G.outerSize(); ++k) {
        for (Eigen::SparseMatrix<double>::InnerIterator it(G,k); it; ++it) {
            out << it.row() << " " << it.col() << " " << it.value() << "\n";
        }
    }
    out.close();
}

int main(){
   std::vector<int> firstPassageTimes;
    // maximum steps for each simulation
    
    bool saveTrajectory = true;
    
    // run simulations
    for (int i = 0; i < numSimulations; ++i) {
        Walker w = {N/2-0.5,N/2+0.5};
        G.setZero();

        int step = 0;
        // run this simulation until the walker reaches the target or until maximum steps
    
        std::ofstream trajectoryFile;
        if (saveTrajectory) {
            trajectoryFile.open("Sim" + std::to_string(i+1) + ".txt");
        }
    
    
    
        while (distanceToTarget(w, t) > a && step < maxSteps) {
            int startX = w.x/1;
            int startY = w.y/1;
            w = jumpCalculate(w);
            if (w.x<0){
                w.x+= (float) N;
            }
            if (w.y<0){
                w.y+=static_cast<float>(N);
            }
            if (w.y>static_cast<float>(N)){
                w.y-=static_cast<float>(N);
            }
            if(w.x>static_cast<float>(N)){
                w.x-=static_cast<float>(N);
            }
            int outX=w.x/1;
            int outY=w.y/1;
           
            if (outX==startX && outY==startY){
                //Do nothing
            } else {
                G.coeffRef(outX,outY) += 1;
            }
            ++step;

            // if switch is on, save the current position
            if (saveTrajectory) {
                trajectoryFile << w.x << " " << w.y << std::endl;
            }
        }

        if (saveTrajectory) {
            trajectoryFile.close();

            writeSparseTxt(G, i+1);
        }

        firstPassageTimes.push_back(step);
        // Print out which simulation we just ran
       // std::cout << "Finished running simulation " << (i+1) << std::endl;
       if ((i+1) % 10 == 0) {
        std::cout << "Finished running simulation " << (i+1) << std::endl;
        }
    }

    double totalFirstPassageTime = 0.0;
    for (int step : firstPassageTimes) {
        totalFirstPassageTime += step;
    }
    double averageFirstPassageTime = totalFirstPassageTime / numSimulations;
    std::cout << "Average first passage time: " << averageFirstPassageTime << std::endl;

    std::ofstream outFile("first_passage_times.txt");
    for (int time : firstPassageTimes) {
        outFile << time << std::endl;
    }
    outFile.close();

return 0;
}

