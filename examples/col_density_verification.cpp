#include <iostream>
#include <fstream>
#include "nuFATE/nuFATE.h"


int main(){

/*
 * flavor: Select a neutrino flavor (1, 2, 3, for electron, muon, and tau, respectively. the number is negative for anti-neutrinos)
 *  gamma: Spectral index of the incoming flux
 *  file: path to file containing the cross sections.
*/
    int flavor = -2;
    double gamma = 2.2;
    bool include_secondaries = false;
    std::string file = "../resources/NuFATECrossSections.h5";

    // ------------------- // ------------------- // ------------------- // ------------------- // ------------------- //
    // normal earth case
    //Initialize an instance of the nuFATE class with these parameters.
    nufate::nuFATE object(flavor, gamma, file, include_secondaries);
    //Result is a struct that stores the solution to the cascade equation.
    nufate::Result result;
    //get_eigs() solves the cascade equation
    result = object.getEigensystem();

    int NumNodes;
    NumNodes = object.getNumNodes();
    std::vector<double> eval = result.eval;
    std::shared_ptr<double> evec = result.evec;
    std::vector<double> ci = result.ci;
    std::vector<double> energy_nodes = result.energy_nodes_;
    std::vector<double> phi_0 = result.phi_0_;
    //Calculate earth column density for a given zenith
    double Na = 6.0221415e23;
    double zenith = 2.2689280276;
    double t;
    t = object.getEarthColumnDensity(zenith) * Na;

    // ------------------- // ------------------- // ------------------- // ------------------- // ------------------- //
    // this is for the shell earth case
    std::vector<double> shellbounds = {0, 0.195, 0.3133, 0.4316, 0.55, 0.6625, 0.775, 0.895, 1};
    std::vector<double> shelldensities = {12.9779, 11.92205, 11.30781, 10.42434, 5.371738, 5.01686, 4.60748, 3.60756};

    nufate::nuFATE object_shell(flavor, gamma, file, include_secondaries, shellbounds, shelldensities);

    nufate::Result result_shell;
    result_shell = object_shell.getEigensystem();

    NumNodes = object_shell.getNumNodes();
    std::vector<double> eval_shell = result_shell.eval;
    std::shared_ptr<double> evec_shell = result_shell.evec;
    std::vector<double> ci_shell = result_shell.ci;
    std::vector<double> energy_nodes_shell = result_shell.energy_nodes_;
    std::vector<double> phi_0_shell = result_shell.phi_0_;
    //Calculate earth column density for a given zenith
    double t_shell;
    t_shell = object_shell.getEarthColumnDensity_shells(zenith) * Na;


    // ------------------- // ------------------- // ------------------- // ------------------- // ------------------- //
    // write out
    std::ofstream outputFile("col_density_comparison.txt"); // Open file for writing (overwrite mode)

    // Check if file is opened successfully
    if (!outputFile) {
        std::cerr << "Error: Unable to open file!" << std::endl;
        return 1; // Return error code
    }
    // Perform the loop and write outputs to the file
    for (double x=1; x>=-1; x=x-0.001){
      std::cout << x << " : " << object.getEarthColumnDensity(std::acos(x)) << " : " << object_shell.getEarthColumnDensity_shells(std::acos(x)) << std::endl;
      outputFile << x << " " << object.getEarthColumnDensity(std::acos(x)) << " " << object_shell.getEarthColumnDensity_shells(std::acos(x)) << std::endl;
    }
    // Close the file
    outputFile.close();

    std::cout << "Written to col_density_comparison.txt" << std::endl;


    // ------------------- // ------------------- // ------------------- // ------------------- // ------------------- //
    // this is for the shell earth case
    std::vector<double> abs;
    std::vector<double> phi_sol;

    //Get Attenuation!
    if(not include_secondaries){
      //Without Secondaries
      for(int i=0; i<NumNodes; i++){
        double sum = 0.;
        for (int j=0; j<NumNodes;j++){
          abs.push_back(ci[j] * exp(t*eval[j]));
          sum+= abs[j] *  *(evec.get()+i*NumNodes+j) / phi_0[i]; // * (std::pow(energy_nodes[i],-2) / std::pow(energy_nodes[i],-gamma));
        }
        phi_sol.push_back(sum);
      }
      //Print Solution
      for(int i =0; i<NumNodes; i++){
        std::cout << phi_sol[i] << std::endl;
      }

    } else{
      //With Secondaries
        int rsize = 2*NumNodes;
        for(int i=0; i<rsize; i++){
          double sum = 0.;
          abs.clear();
          for (int j=0; j<rsize;j++){
            abs.push_back(ci[j] * exp(-t*eval[j]));
            sum+= (abs[j] * *(result.evec.get()+i*rsize+j)) / phi_0[i];
          }
        phi_sol.push_back(sum);
        }

        //Print Solution
        std::cout << "Solution Including secondaries= " << std::endl;
        for(int i =0; i<NumNodes; i++){
          std::cout << phi_sol[i] << std::endl;
        }
    }

    return 0;
}
