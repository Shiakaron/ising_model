#include <iostream>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <thread>
#include <string>
#include <sstream>
#include <future>
//testing these
#include <thread>
#include <mutex>

using namespace std;

struct quantity{
    double value;
    double error;
};

extern int L;          //size of array
extern int N;                //number of sites
extern int **F;                    //spin array pointer (2dim)
extern unsigned int seed;
extern int *S;          //wolf spin array (n dimensional)
extern int *F_2;         //metropolis spin array (n dimensional)

extern int nCycles;       //number of cycles of Monte Carlo primary
extern int mCycles;          //number of cycles of Monte Carlo secondary
extern int dataPoints;     //total data points

extern int n_b;               //jackknife error analysis. n_b is number of bins

extern double iniT;        //initial temperature
extern double finT;        //final temperature
extern int numT;           //data points (for linear spacing)
extern double *T;                 //array of temperatures

extern double tempM;        //temporary magnetisation

extern double Next2Nearest; // next to nearest neighbour interaction
extern double B; // external B field
extern int dim;     // number of dimensions for generalised metropolis and wolff
extern int V;       // Volume of n dimensional lattice

//GENERAL USE
double* linspace(double a, double b, int N);
double sumArray(double arr[], int siz);
double averageArray(double arr[], int siz);
void printArray();
void printArraytoFile(ofstream &file);
void printArraytoFile2(ofstream &file);
void printArray3D(int arr[]);
void printArray3D_2(int arr[], const string &filename);
void printArraytoFile3D(ofstream &file);
void printArraytoFile3D_2(const string &filename, double T);
bool file_exists(const string& name);
void filename_rename_if_exists(string& filename);
void print_T_to_console();

//INITIALISATIONS
void initiStateCold();
void initiStateHot();

//LATTICE ESTIMATIONS
double denergy(int i, int j);

//double flipFunction(double T, int cycles);
void flipFunction(double T, int cycles);
void flipFunction2(double T, int cycles);
double computeEnergy();
double computeMagnetisation();
double computeSusceptibility(double arr[], int siz, double T);
double chiArray(double arrayM[], int binSize, double T);
double calcSigma(double arr[], int siz);
double jackknife(int n_b, double arr[]) ;
double* jackknife2(int n_b, double arr[]);
double jackknife_new(double data[], int n_b);

//ANALYSIS LINES
void jackknifeErrorAnalysis();
void jackknifeErrorAnalysis2();
void MagnetisationEnergyData();
void MagnetisationEnergyDataNew();
void ThermalisationInvestigation();
void simulationOfSpins(double T, int sCycles, int nConfigurations);
void simulationOfSpins2(double T, int sCycles, int nConfigurations);
void ThermalisationInvestigationNew();
void susceptibilityData();
void susceptibilityData2(int L_arg);
quantity susceptibilityData3(double T_mt);
void testMetropolisND();
void testMetropolisND2();
void testMetropolisND2_1(double T, int cycles, int configs);
void produceConfigurations3D(int L_arg, int cycles, int configs);
void produceConfigurations3D_2(double T, int cycles, int configs);
quantity magnetisationData_mt(double T_mt);

//AUTOCORRELATION
void autocorrelationInvestigationNew();
void autocorrelationInvestigation();
double rho_calculation(double arr[], int t, int t_max);
double* autocorrelation(double arr[], int t_max);
double jackknife_new(double data[], int n_b);

//WOLFF
quantity wolff(double T);
void init_hot_wolff(int arr[]);
void print_2D_array(int arr[]);
double average_grain_size(int *S);
double bootstrap_error(double data[], int array_size, int bin_number);

//METROPOLIS_N_DIMENSIONS
void initStateColdND();
void initiStateHotND();
double denergyND(int index);
void flipFunctionND(double T, int cycles);
double computeMagnetisationND();
//for multithreading
void initiStateHotND_mt(int arr[]);
double computeMagnetisationND_mt(int arr[]);
double denergyND_periodic_mt(int arr[], int index);
void flipFunctionND_mt(int arr[], double T_mt, double& tempM_mt,  int cycles);
