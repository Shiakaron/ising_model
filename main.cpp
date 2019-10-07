#include "header.h"

int dim;
int L;
int N;
int *spins;
map<int,vector<int>> mapOfNearest;
map<int,vector<int>> mapOfNext2Nearest;
int M;

double n2n = 0;
double H = 0;

unsigned int seed = (unsigned)time(0); // (unsigned)time(0);
double *T;

void magnetisation_vs_time_data()
{
    cout << "Running for magnetisation of evolving system" << endl;
    L = 5;
    dim = 2;
    int thermalisationCycles = 0;
    int dataPoints = 5000;
    int spacingCycles = 1;
    double Temp = 5.0;

    // initialise vector
    vector<double> magn;

    //initialise spins array and neighbours maps
    initialise_system_and_maps();
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);

    // open file
    ofstream myfile;
    string folder = ".\\data\\magn_vs_time";
    string filename = "magn_vs_time_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(Temp)+".txt";
    string path = folder+"\\"+filename;
    filename_rename_if_exists(filename, folder);
    myfile.open(path);
    if (!myfile.is_open()) {
        throw path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // begin
    initialise_spins_hot();
    compute_magnetisation();
    metropolis_function(Temp,thermalisationCycles);
    magn.push_back(fabs(M));
    for (int i=1; i<dataPoints; i++) {
        metropolis_function(Temp,spacingCycles);
        magn.push_back(fabs(M));
    }

    // write data on opened file
    myfile << Temp << endl;
    for(int i = 0; i < dataPoints; i++) {
      myfile << magn[i] << endl;
    }

    // wrap up
    delete[] spins;
    if (myfile.is_open()){
       myfile.close();
    }
}

void magnetisation_vs_temp_data()
{
    cout << "Running for magnetisation at different temperatures data" << endl;
    L = 10;
    dim = 2;
    int thermalisationCycles = 5000;
    int spacingCycles = 50;
    int dataPoints = 200;      //total data points
    double iniT = 1.0; double finT = 5.0; int numT = 41;
    T = linspace(iniT, finT, numT);

    // initialise vectors, pointers and time variable
    vector<double> magn;
    vector<double> err_magn;
    double *arrayM;
    double *bootstrap_values;
    clock_t tStartTemp;

    //initialise spins array and neighbours maps
    initialise_system_and_maps();
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);

    // open file
    ofstream myfile;
    string folder = ".\\data\\magn_data";
    string filename = "magn_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    string path = folder+"\\"+filename;
    filename_rename_if_exists(filename, folder);
    myfile.open(path);
    if (!myfile.is_open()) {
        throw path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        cout << "T = " << T[i] << endl;
        tStartTemp = clock();
        arrayM = new double[dataPoints];
        // begin
        initialise_spins_hot();
        compute_magnetisation();
        metropolis_function(T[i],thermalisationCycles);
        arrayM[0] = fabs((double)M/N);
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            arrayM[j] = fabs((double)M/N);
        }
        // average and error
        bootstrap_values = bootstrap_error(arrayM, dataPoints, 128, true);
        magn.push_back(bootstrap_values[0]);
        err_magn.push_back(bootstrap_values[1]);
        // wrap up
        delete[] arrayM;
        delete[] bootstrap_values;
        cout << "Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

    // write data on opened file
   for(int i = 0; i < numT; i++) {
      myfile << T[i] << " " << magn[i] << " " << err_magn[i] << endl;
   }

   // wrap up
   delete[] spins;
   delete[] T;
   if (myfile.is_open()){
       myfile.close();
   }
}

int main()
{
    clock_t tStart = clock();
    srand(seed);
    try {
        // magnetisation_vs_time_data();
        magnetisation_vs_temp_data();
        // investigating_autocorrelation();
    }
    catch (string path){
        cerr << "ERROR: FILE NOT OPENED. EXITING PROGRAM." << endl;
        cerr << "Attempted file path: " << path << endl;
        return 1;
    }
    cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}
