#include <iostream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <map>
#include <vector>

using namespace std;

// struct observable{
//     double value;
//     double error;
// };

extern int dim;
extern int L;
extern int N;
extern int *spins;
extern map<int,vector<int>> mapOfNearest;
extern map<int,vector<int>> mapOfNext2Nearest;
extern int M;
extern double n2n;
extern double H;
extern unsigned int seed;
extern double *T;
extern int dataPoints;

//MAIN
void magnetisation_vs_temp_data();
void magnetisation_vs_time_data();

//METROPOLIS
void initialise_system_and_maps();
void initialise_spins_cold();
void initialise_spins_hot();
void initialise_nearest_periodic();
void initialise_next2nearest_periodic();
double compute_denergy();
void compute_magnetisation();
void metropolis_function(double T, int cycles);

//GENERAL
double* linspace(double a, double b, int N);
double sumArray(double arr[], int siz);
double averageArray(double arr[], int siz);
bool file_exists(const string& p_file);
void filename_rename_if_exists(string& filename, string& folder);
void print_T(int siz);
void print_mapOfNearest();
void print_mapOfNext2Nearest();
void print_all_parameters(int thermalisationCycles, int dataPoints, int spacingCycles, int numT, double Temp);
void print_spins();
void print_array(double arr[], int siz);
void print_spins_2D();

//BOOTSTRAP ERROR ANALYSIS
double* bootstrap_error(double data[], int array_size, int bin_number, bool decorelate);

//AUTOCORRELATION
double* autocorrelation(double arr[], int t_max);
double rho_calculation(double arr[], int t, int t_max);
