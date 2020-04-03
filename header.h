#include <iostream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <limits>

using namespace std;

struct observable{
    double value;
    double error;
};

extern int dim;
extern int L;
extern int N;
extern double N_links;
extern int *spins;
extern map<int,vector<int>> mapOfNearest;
extern map<int,vector<int>> mapOfNext2Nearest;
extern int M;
extern double n2n;
extern double H;
extern unsigned int seed;
extern double *T;
extern int dataPoints;
extern double E;

//MAIN
int initial_menu();
void magnetisation_vs_time_data();
void magnetisation_vs_time_data_bulk();
void magnetisation_vs_temp_data();
void autocorrelation_initial_investigation();
void autocorrelation_peak_investigation();
void energy_first_check();
void energy_second_check();
void energy_vs_time_data();
void energy_vs_temp_data();
void heat_capacity_data();
void heat_capacity_peak_data();
void external_field_investigation();
void magnetic_susceptibility_peak_data();
void generate_configurations_for_gif();
void generate_configuration_for_figure();
void next_to_nearest_investigation();
void wolff_cluster_size_vs_temp_data();

//METROPOLIS
void initialise_system_and_maps();
void initialise_spins_auto(double Temp);
void initialise_spins_cold();
void initialise_spins_hot();
void initialise_nearest_periodic();
void initialise_next2nearest_periodic();
double compute_denergy();
void compute_magnetisation();
void metropolis_function(double Temp, int cycles);
void compute_energy();
double energy_per_link();
double energy_per_site();
double* get_heat_capacity(double arr[], int siz, double Temp);
double* get_magnetic_susceptibility(double arr[], int siz, double Temp);

//WOLFF
double* wolff_function(double Temp, int cycles);
int flip_cluster(int flags[], int seed_spin);

//GENERAL
double* linspace(double a, double b, int N);
bool file_exists(const string& p_file);
void filename_rename_if_exists(string& filename, string& folder);
double sumArray(double arr[], int siz);
double averageArray(double arr[], int siz);
observable compute_average_and_sigma(double arr[], int siz);
void print_T(int siz);
void print_mapOfNearest(int max);
void print_mapOfNext2Nearest(int max);
void print_all_parameters(int thermalisationCycles, int dataPoints, int spacingCycles, int numT, double Temp);
void print_spins();
void print_array(double arr[], int siz);
void print_spins_2D();
int user_integer_input(int min, int max);

//BOOTSTRAP ERROR ANALYSIS
double* bootstrap_error(double data[], int array_size, int bin_number, bool decorelate);

//AUTOCORRELATION
double* autocorrelation(double arr[], int t_max);
double rho_calculation(double arr[], int t, int t_max);
