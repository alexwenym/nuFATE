#include <vector>
#include <string>
#include <math.h>
#include <cmath>
#include <memory>
#include <iostream>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "H5Gpublic.h"
#include "H5Fpublic.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_integration.h>

#ifndef NUFATE_H
#define NUFATE_H

namespace nufate{

/// \struct Square_matrix_double
/// \brief Very simple matrix struct used for pybinding.
/// vec_ is the front address of a chank of double data,
/// which will be shaped in NxN 2dim matrix in python.
struct Square_matrix_double {
 public :
 unsigned int n_;
 std::shared_ptr<double> vec_;
};

///\class Result
///\brief Simple class to hold the results.
class Result {
  public :
   std::vector<double> eval;          // eigen values
   std::shared_ptr<double> evec;      // pointer to top address of eigen vectors
   std::vector<double> ci;            // coefficients
   std::vector<double> energy_nodes_; // energy nodes in GeV
   std::vector<double> phi_0_;        // initial flux * E^(pedestal_gamma)

   //
   // functions below are not used by main nuFATE program,
   // but may be used for pybinding.
   //

   /// \brief getter for eigenvalues
   /// @return 1D vector
   std::vector<double> get_eigenvalues() { return eval; }

   /// \brief getter for coefficients
   /// @return 1D vector
   std::vector<double> get_coefficients()   { return ci; }

   /// \brief getter for energy nodes
   /// @return 1D vector
   std::vector<double> get_energy_nodes() { return energy_nodes_; }

   /// \brief getter for initial flux * E^(pedestal_gamma)
   /// @return 1D vector
   std::vector<double> get_phi_0() { return phi_0_; }

   /// \brief converter from pointer of double to Square_matrix_double object
   /// @return Square_matrix_double
   Square_matrix_double get_eigenvec_matrix() { 
      unsigned int n = energy_nodes_.size();
      Square_matrix_double smatrix;
      // evec is n x n data
      smatrix.n_ = n;
      smatrix.vec_ = evec;   
      return smatrix;
   }

   /// \brief debug print function
   void Print(unsigned int target_index = 0) {
      unsigned int dim = energy_nodes_.size();

      std::cout << "***** Result print *****"  << std::endl;
      std::cout << "eigenvectors at row " << target_index << " ---"  << std::endl;
      for (unsigned int i=target_index; i<target_index+1; ++i) {
         for (unsigned int j=0; j<dim; ++j) {
             std::cout << *(evec.get() + i*dim + j) << " " ;
         }
         std::cout << std::endl;
      }

      std::cout << std::endl;
      std::cout << "eigenvalue at i = " << target_index << " ---"  << std::endl;
      std::cout << eval[target_index] << std::endl;
 
      std::cout << std::endl;
      std::cout << "coefficient at i = " << target_index << " ---"  << std::endl;
      std::cout << ci[target_index] << std::endl;

      std::cout << std::endl;
      std::cout << "energy_node at i = " << target_index << " ---"  << std::endl;
      std::cout << energy_nodes_[target_index] << std::endl;
   }
};


///\class nuFATE
///\brief nuFATE main class
class nuFATE {
  private:
    const double GF = 1.16e-5;
    const double hbarc = 1.97e-14;
    const double GW = 2.085;
    const double MW = 80.385e0;
    const double mmu = 0.106e0;
    const double me = 511.e-6;
    const double pi = M_PI;
    const double MZ = 91.18;
    const double s2t = 0.23;
    const double gL =  s2t-0.5;
    const double REarth = 6371.;
    const double gR = s2t;
  private:
    int newflavor_;
    double newgamma_;
    std::string energy_filename; 
    std::string newh5_filename_;
  public:
    /// \brief Constructor
    /// @param flavor position ID of the system (1:NuE, -1:NuEBar, 2:NuMu, -2:NuMuBar, 3:NuTau, -3:NuTauBar)
    /// @param gamma spectral index of input flux.
    /// @param h5_filename name of hdf5 file containing the cross sections.
    /// @param include_secondaries if true secondaries are added to the flux propagation.
    nuFATE(int flv, double gamma, std::string h5_filename, bool include_secondaries);
    /// \brief Constructor
    /// @param flavor position ID of the system (1:NuE, -1:NuEBar, 2:NuMu, -2:NuMuBar, 3:NuTau, -3:NuTauBar)
    /// @param gamma spectral index of input flux.
    /// @param h5_filename name of hdf5 file containing the cross sections.
    /// @param include_secondaries if true secondaries are added to the flux propagation.
    nuFATE(int flv, double gamma, std::string h5_filename, bool include_secondaries, const std::vector<double>& radialbounds, const std::vector<double>& radialdensities);
    // /// \brief Constructor
    // /// @param flavor position ID of the system (1:NuE, -1:NuEBar, 2:NuMu, -2:NuMuBar, 3:NuTau, -3:NuTauBar)
    // /// @param energy_filename file containing the desired energy spectrum
    // /// @param h5_filename name of hdf5 file containing the cross sections.
    // /// @param include_secondaries if true secondaries are added to the flux propagation.
    // nuFATE(int flv, std::string energy_filename, std::string h5_filename, bool include_secondaries);
    /// \brief Constructor
    /// @param flavor position ID of the system (1:NuE, -1:NuEBar, 2:NuMu, -2:NuMuBar, 3:NuTau, -3:NuTauBar)
    /// @param gamma spectral index of input flux.
    /// @param energy_nodes energy nodes in GeV.
    /// @param sigma_array total cross section(CC+NC) at each energy node in cm^2.
    /// @param dsigma_dE square array of differential cross section (NC only) in cm^2/GeV.
    /// @param include_secondaries if true secondaries are added to the flux propagation.
    nuFATE(int flv, double gamma, std::vector<double> energy_nodes, std::vector<double> sigma_array, std::vector<std::vector<double>> dsigma_dE, bool include_secondaries);
    /// \brief Constructor
    /// @param flavor position ID of the system (1:NuE, -1:NuEBar, 2:NuMu, -2:NuMuBar, 3:NuTau, -3:NuTauBar)
    /// @param gamma spectral index of input flux.
    /// @param cc_sigma_file file path for charge current(CC) total cross section in input format of nuSQuIDS
    /// @param nc_sigma_file file path for neutral current(NC) total cross section in input format of nuSQuIDS
    /// @param nc_dsigma_file file path to NC dSigma/dE differential cross section in input format of nuSQuIDS
    /// @param include_secondaries if true secondaries are added to the flux propagation.
    nuFATE(int flv, double gamma, std::string cc_sigma_file, std::string nc_sigma_file, std::string nc_dsigma_file, bool include_secondaries);
 
    /// \brief Eigensystem calculator
    Result getEigensystem();
    /// \brief Function to get Earth column density
    /// @param theta Zenith angle in radians.
    double getEarthColumnDensity(double theta);
    /// \brief Function to get Earth column density given a shell model with specified list of bounds and const. densities
    double getEarthColumnDensity_shells(double theta);
    /// \brief Function to get flavor
    int getFlavor() const;
    /// \brief Function to get input spectral index
    double getGamma() const;
    /// \brief Function to get input filename
    std::string getFilename() const;
    /// \brief Function to get number of energy nodes
    double getNumNodes() const;
    /// \brief Function to toggle secondaries
    void setAddSecondaries(bool opt) { add_secondary_term_ = opt;}

    /// \brief Function to get energy nodes
    std::vector<double> getEnergyNodes() { return energy_nodes_; }
    /// \brief Function to get total crosssections (CC + NC)
    std::vector<double> getTotalCrossSections() { return sigma_array_; }
    /// \brief Function to get dsigma_dE differential crossection (for neutral-current only)
    Square_matrix_double getNCDifferentialCrossSections() {
       Square_matrix_double smatrix;
       smatrix.n_ = NumNodes_;
       smatrix.vec_ = dxs_array_;
       return smatrix; 
    }

    /// \brief Calculate Relative Attenuation 
    /// @param total number of targets, X[g/cm^2]*Na(Avogadro number) for example
    /// @return attenuation factor (arrival_flux / initial_flux), 1D array in energy bins
    std::vector<double> getRelativeAttenuation(double number_of_targets);

  protected:
    void Init(double gamma, const std::vector<double> &enodes, const std::vector<double> &sigmas, const std::vector<std::vector<double> > &dsigma);
    void AddAdditionalTerms();
    void LoadCrossSectionFromHDF5();
    void SetCrossSectionsFromInput(const std::vector<double> &sigma, const std::vector<std::vector<double>> &dsigma_dE);
    void SetInitialFlux();
    void set_glashow_total();
    void set_glashow_partial();
    void set_RHS_matrices(std::shared_ptr<double> RMatrix_, std::shared_ptr<double> dxs_array);
    static double rho_earth(double theta, void * p); //returns density for a given angle in radians along x distance
    double readDoubleAttribute(hid_t, std::string) const;
    unsigned int readUIntAttribute(hid_t, std::string) const;
    std::vector<double> logspace(double min,double max,unsigned int samples) const;
    double chord(double perpdist, double radius); 
    void set_shellEarth(const std::vector<double> &radialbounds, const std::vector<double> &radialdensities);

  private:
    void AllocateMemoryForMembers(unsigned int num_nodes);
    void SetEnergyBinWidths();
    std::string grpdiff_;
    std::string grptot_;
    hid_t file_id_;
    hid_t root_id_;
    unsigned int NumNodes_;
    unsigned int rsize_;
    double Emax_;
    double Emin_;
    int dxsdim_[2];
    std::vector<double> energy_nodes_;
    std::vector<double> sigma_array_;
    std::vector<double> sigma_array_orig_;
    std::vector<double> DeltaE_;
    std::vector<double> phi_0_;
    std::vector<double> glashow_total_;
    std::vector<double> sig3_array_;

    std::vector<double> shell_radialbounds_;
    std::vector<double> shell_angularbounds_;
    std::vector<double> shell_radialdensities_;

    std::shared_ptr<double> glashow_partial_;
    std::shared_ptr<double> dxs_array_;
    std::shared_ptr<double> tau_array_;
    std::shared_ptr<double> sec_array_;
    std::shared_ptr<double> regen_array_;
    std::shared_ptr<double> RHSMatrix_;
    std::shared_ptr<double> RHSMatrix1_;
    std::shared_ptr<double> RHSMatrix2_;
    std::shared_ptr<double> RHSMatrix3_;
    std::shared_ptr<double> RHSMatrix4_;
    std::shared_ptr<double> RHregen_;
    std::shared_ptr<double> Enuin_;
    std::shared_ptr<double> Enu_;
    std::shared_ptr<double> selectron_;
    std::shared_ptr<double> den_;
    std::shared_ptr<double> t1_;
    std::shared_ptr<double> t2_;
    std::shared_ptr<double> t3_;
  private:
    bool include_secondaries_ = false;
    bool memory_allocated_ = false;
    bool initial_flux_set_ = false;
    bool total_cross_section_set_ = false;
    bool differential_cross_section_set_ = false;
    bool RHS_set_ = false;
    bool add_tau_regeneration_ = true;
    bool add_glashow_term_= true;
    bool add_secondary_term_ = true;

    bool use_shell_earth_ = false; 
};

}

#endif
