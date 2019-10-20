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

/*
To start with I will analyse the evolution of the average magnetisation of the system.
For different Lattice sizes and temperatures I will output different data files.
PLOT 1: For a fixed lattice size plot the evolution of the system at different
temperatures, both below and above Tc.
PLOT 2: I will choose a specific lattice size and temperature (below Tc)
to demonstrate the randomness of the simulation and difference of hot and cold start
PLOT 3: For a fixed temperature (below Tc) plot the evolution of the system at
different lattice sizes. Used average magnetisation for clarity in plots.
*/

void magnetisation_vs_time_data()
{
    L = 40;
    dim = 2;
    int thermalisationCycles = 0;
    int dataPoints = 5000;
    int spacingCycles = 1;
    double Temp = 1.5;

    //initialise spins array and neighbours maps
    initialise_system_and_maps();
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);
    cout << "Running for magnetisation of evolving system" << endl;

    // open file
    ofstream myfile;
    string folder = ".\\data\\magn_vs_time";
    string filename = "magn_vs_time_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(Temp)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: magnetisation_vs_time_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;

    // begin
    cout << "Starting computations" << endl;
    myfile << Temp << endl;
    // initialise_spins_hot();
    initialise_spins_cold();
    compute_magnetisation();
    metropolis_function(Temp,thermalisationCycles);
    myfile << (fabs(M)) << endl;
    for (int i=1; i<dataPoints; i++) {
        metropolis_function(Temp,spacingCycles);
        myfile << (fabs(M)) << endl;
    }

    // wrap up
    delete[] spins;
    if (myfile.is_open()){
       myfile.close();
    }
}

void magnetisation_vs_time_data_bulk()
{
    dim = 2;
    int L_list[5] = {5,10,20,40,80};
    int thermalisationCycles = 0;
    int dataPoints = 5000;
    int spacingCycles = 1;
    double iniT = 1.0; double finT = 5.0; int numT = 9;
    T = linspace(iniT, finT, numT);

    for (int i=0; i<numT; i++) {
        for (int j=0; j<5; j++) {
            //set L
            L = L_list[j];
            //initialise vector
            vector<double> magn;

            //initialise spins array and neighbours maps
            initialise_system_and_maps();
            print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, T[i]);
            cout << "Running for magnetisation of evolving systems" << endl;
            //open file
            ofstream myfile;
            string folder = ".\\data\\magn_vs_time";
            string filename = "magn_vs_time_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(T[i])+".txt";
            filename_rename_if_exists(filename, folder);
            string path = folder+"\\"+filename;
            myfile.open(path);
            if (!myfile.is_open()) {
                throw "Func: magnetisation_vs_time_data_bulk(). File not opened with path: "+path;
            }
            cout << "Writing in file with path: " << path << endl;
            cout << "Starting computations" << endl;
            //begin
            initialise_spins_hot();
            compute_magnetisation();
            metropolis_function(T[i],thermalisationCycles);
            magn.push_back(fabs(M));
            for (int k=1; k<dataPoints; k++) {
                metropolis_function(T[i],spacingCycles);
                magn.push_back(fabs(M));
            }

            // write data on opened file
            myfile << T[i] << endl;
            for(int k=0; k<dataPoints; k++) {
              myfile << magn[k] << endl;
            }

            // wrap up
            delete[] spins;
            if (myfile.is_open()){
               myfile.close();
            }
        }
    }
}

void magnetisation_vs_temp_data()
{
    L = 64;
    dim = 2;
    int thermalisationCycles = 10000;
    int spacingCycles = 50;
    int dataPoints = 2000;      //total data points
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
    cout << "Running for magnetisation at different temperatures data" << endl;
    // open file
    ofstream myfile;
    string folder = ".\\data\\magn_data";
    string filename = "magn_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: magnetisation_vs_time_data(). File not opened with path: "+path;
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
        magnetisation_vs_time_data();
        // magnetisation_vs_time_data_bulk();
        // magnetisation_vs_temp_data();
    }
    catch(string er) {
        cerr << "ERROR. EXITING PROGRAM." << endl;
        cout << er << endl;
        return 1;
    }
    cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}
