#include "header.h"

int dim; // Dimensions
int L; // Lattize size
int N; // Total number of sites
double N_links; // Number of edges
int *spins; // poitner to spin array
map<int,vector<int>> mapOfNearest; // map storing all nearest neighbours for each site
map<int,vector<int>> mapOfNext2Nearest; // map storing all next-to-nearest neigbours for each site
int M; // magnetisation of the system
double E; // energy of the system

double n2n = 0; // next-to-neareset neighbour interaction
double H = 0; // external field strength

unsigned int seed = (unsigned)time(0); // (unsigned)time(0); Randomized seed according to time on operating system.
double *T; // pointer to array of temperatures over which our analysis will iterate

int main(int argc, char** argv)
{
    clock_t tStart = clock();
    srand(seed);
    cout << "Ising Model, made for Part II Physics Computing Project" << endl;
    cout << "C++ was used to produce data files ('low' level -> faster computation) which were later analysed and plotted using Python.\n";

    // user friendly function to run the code. The user will input an integer to choose the desired analysis line. This is enclosed in a try-catch block in case there is a run-time error.
    int user_choice = initial_menu();
    try {
        switch(user_choice) {
            case 1: magnetisation_vs_time_data();
            break;
            case 2: magnetisation_vs_temp_data();
            break;
            case 3: autocorrelation_initial_investigation();
            break;
            case 4: autocorrelation_peak_investigation();
            break;
            case 5: energy_vs_time_data();
            break;
            case 6: energy_vs_temp_data();
            break;
            case 7: heat_capacity_data();
            break;
            case 8: heat_capacity_peak_data();
            break;
            case 9: external_field_investigation();
            break;
            case 10: magnetic_susceptibility_peak_data();
            break;
            case 11: generate_configurations_for_gif();
            break;
            case 12: generate_configuration_for_figure();
            break;
            case 13: next_to_nearest_investigation();
            break;
            case 14: wolff_cluster_size_vs_temp_data();
            break;
            case 15: wolff_autocorrelation_investigation();
            break;
        }
        cout << "\nOperation complete.\n";
    }
    catch(string er) {
        cerr << "\nERROR. EXITING PROGRAM." << endl;
        cout << er << endl;
        return 1;
    }
    cout << "\nProgram ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}

int initial_menu()
// User friendly function to choose analysis line
{
    cout << "\nPlease select an analysis by entering an integer:\n";
    cout << "1 for Magnetisation vs Time.\n";
    cout << "2 for Magnetisation vs Temperature.\n";
    cout << "3 for Tau_e vs Temperature initial investigation.\n";
    cout << "4 for Tau_e vs Temperature close to critical temperature.\n";
    cout << "5 for Energy vs Time.\n";
    cout << "6 for Energy vs Temperature.\n";
    cout << "7 for Heat capacity vs Temperature.\n";
    cout << "8 for Heat capacity vs Temperature around critical point.\n";
    cout << "9 for Magnetisation vs External Field.\n";
    cout << "10 for Magnetic Susceptibility vs Temperature around critical point.\n";
    cout << "11 to generate configurations for a GIF.\n";
    cout << "12 to generate configuration for a Figure.\n";
    cout << "13 for Magnetisation and Energy vs Temperature with Next to Nearest interactions.\n";
    cout << "14 for Cluster size vs Temperaure using Wolff's algorithm.\n";
    cout << "15 for Wolff's algorithm autocorrelation investigation.\n";
    cout << "0 to Exit the program.\n";
    int choice = user_integer_input(0,15);
    return choice;
}

/*
TESTING THE PROGRAM'S CORE FUNCTIONS
To start with analysing the evolution of the average magnetisation of the system. I will do this with different Lattice sizes and a range of temperatures. Following that I will test how the average magnetisation of an evolved system changes with temperature. I will do this for different lattice sizes to show finite size effects.

PLOT 1:
For a fixed lattice size plot the evolution of the system at different temperatures, below, at and above Tc.

PLOT 2:
I will choose a specific lattice size and temperature (below Tc) to demonstrate the randomness of the simulation and difference between hot and cold start

PLOT 3:
For a fixed temperature (below Tc) plot the evolution of the system at different lattice sizes. Use average magnetisation for clarity in plots.

PLOT 4:
<|magnetisaiton|> vs temperature for different lattice sizes. Shows finite size effects and how greater lattice sizes are closer to the real solution.
*/

void magnetisation_vs_time_data()
{
    cout << "Running for magnetisation of evolving system" << endl;
    // set parameters
    dim = 2;
    L = 40;
    int thermalisationCycles = 0; // Monte Carlo Sweeps (MCS) to thermalise the system
    int spacingCycles = 1; // MCS between data points
    int dataPoints = 5000;
    double Temp = 1.5;
    // print out all parameters to the user to check
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);
    // check if the parameters are "OK". give option to change them
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,4]:\n";
        dim = user_integer_input(2,4);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Next to nearest interaction [0,200]/100:\n";
        n2n = user_integer_input(0,200);
        cout << "External field [0,200]/100:\n";
        H = user_integer_input(0,200);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Temperature [1,300]/10:\n";
        Temp = double(user_integer_input(1,300))/10;
    }
    // choose how to initialise the system
    cout << "How to initialise the system: Enter 1 for Hot Start, 0 for Cold Start\n";
    int start = user_integer_input(0,1);

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // initialise: hot or cold
    if (start == 1) {initialise_spins_hot();}
    else {initialise_spins_cold();}

    // open file to write data
    ofstream myfile;
    string folder = ".\\data\\magn_vs_time"; // folder path
    string filename = "magn_vs_time_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(Temp)+".txt"; // file name
    filename_rename_if_exists(filename, folder); // function to check if a file with the same name already exists. if it does add (number) to end of file name
    string path = folder+"\\"+filename; // full path for file
    myfile.open(path);
    if (!myfile.is_open()) {
        // if the file is not opened throw an error for the user to check the specified path
        throw "Func: magnetisation_vs_time_data(). File not opened with path: " + path + "\nPlease fix path";
    }
    cout << "Writing in file with path: " << path << endl;

    // begin computing
    cout << "Starting computations" << endl;
    compute_magnetisation(); // start by copmuting M. this will be updated automatically as the system evolves
    metropolis_function(Temp,thermalisationCycles); // call the metropolis function to evolve the system at temperature Temp and to thermalise the system.
    myfile << fabs((double)M/N) << endl; // output the <|magnetisation|>
    for (int i=1; i<dataPoints; i++) {
        metropolis_function(Temp,spacingCycles); // record data every spacingCycles
        myfile << fabs((double)M/N) << endl;
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
    dim = 2;
    L = 40;
    int thermalisationCycles = 1000;
    int spacingCycles = 50;
    int dataPoints = 2000;      //total data points
    double iniT = 1.0; double finT = 5.0; int numT = 41;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,4]:\n";
        dim = user_integer_input(2,4);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Next to nearest interaction [0,200]/100:\n";
        n2n = user_integer_input(0,200);
        cout << "External field [0,200]/100:\n";
        H = user_integer_input(0,200);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // initialise pointers and time variable
    double *arrayM;
    double *bootstrap_values;
    clock_t tStartTemp;
    arrayM = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\magn_data";
    string filename = "magn_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: magnetisation_vs_temp_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // begin
        initialise_spins_auto(T[i]);
        compute_magnetisation();
        metropolis_function(T[i],thermalisationCycles);
        arrayM[0] = fabs((double)M/N);
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            arrayM[j] = fabs((double)M/N);
        }
        // average and error
        bootstrap_values = bootstrap_error(arrayM, dataPoints, 128, true);
        // write data in myfile
        myfile << T[i] << " " << bootstrap_values[0] << " " << bootstrap_values[1] << endl;
        cout << "T = " << T[i] << ", m = " << bootstrap_values[0] << " +- " <<  bootstrap_values[1] << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

   // wrap up
   delete[] spins;
   delete[] T;
   delete[] arrayM;
   delete[] bootstrap_values;
   if (myfile.is_open()){
       myfile.close();
   }
}

/*
AUTOCORRELATION
Firstly, autocorrelation_initial_investigation() will investigate the autcorrelation function and it's behaviour over a wide range of temperatures and over a "long time". A single value will be recorded at each temperature where the autocorr drops by 1/e, tau_e.

Secondly I will focus my attention around the peak, T_c ~ 2.27 K, where I will attempt to get an estimate of the time lag tau_e by computing the value multiple times.

In the paper M. P. Nightingale and H. W. J. Blote, Phys. Rev. Lett. 76 (1996) they found: "from a finite-size scaling analysis of these autocorrelation times, the dynamic critical exponent z is determined as z = 2.1665 (12)" where "tau_e is similar to L^z at the incipient critical point."

PLOT 1:
Tau_e vs temperature for range 1-5 Kelvin. I will collect data for a range of Lattice sizes. This will give rough idea of what to expect

PLOT 2,3:
For each L and T compute multiple tau_e's around critical point to get average (up to the point where computational time is too much of a cost). Fit a gaussian on each L to get the peak values. Then plot tau_e_peak vs Lattice size and fit a*L^z to get a relationship.
*/

void autocorrelation_initial_investigation()
{
    cout << "Running autocorrelation initial investigation" << endl;
    L = 16;
    dim = 2;
    int thermalisationCycles = 15000; // enough to ensure thermalisation
    int dataPoints = 5000; // "long time"
    int spacingCycles = 1;
    double iniT = 1.5; double finT = 4.5; int numT = 31;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,3]:\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(1000,20000);
        cout << "Initial Temperature [1,49]/10:\n";
        int iniT_int = user_integer_input(1,49);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",50]/10:\n";
        int finT_int = user_integer_input(iniT_int,50);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    //initialise vectors, pointers and time variable
    double *arrayM;
    double *bootstrap_values;
    double *autocorr;
    int tau_e;
    clock_t tStartTemp;
    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open myfile -> record: T , tau_e
    ofstream myfile;
    string folder = ".\\data\\autocorrelation_data\\initial_investigation";
    string filename = "autocorr_times_"+to_string(L)+"_L.txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: autocorrelation_initial_investigation(). File not opened with path: " + path + "\nPlease fix path";
    }
    cout << "Writing in file with path: " << path << endl;

    // myfile2 -> bulk data, for each L and T record: magn, autocorr
    ofstream myfile2;
    string folder2 = ".\\data\\autocorrelation_data\\initial_investigation\\bulkdata";
    string filename2;
    string path2;

    // begin
    arrayM = new double[dataPoints];
    cout << "Starting computations" << endl;
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // open myfile2
        filename2 = "autocorr_times_"+to_string(L)+"_L_"+to_string(T[i])+".txt";
        filename_rename_if_exists(filename2, folder2);
        path2 = folder2+"\\"+filename2;
        myfile2.open(path2);
        if (!myfile2.is_open()) {
            throw "Func: autocorrelation_initial_investigation(). File not opened with path: " + path2 + "\nPlease fix path";
        }

        // magnetisation array
        initialise_spins_cold();
        compute_magnetisation();
        metropolis_function(T[i],thermalisationCycles);
        arrayM[0] = fabs((double)M/N);
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            arrayM[j] = fabs((double)M/N);
        }

        // autocorrelation array
        autocorr = autocorrelation(arrayM, dataPoints);
        tau_e = dataPoints;
        for (int j = 0; j < dataPoints; j++){
           if (autocorr[j] <= exp(-1)){
               tau_e = j;
               break;
           }
        }

        // write in myfile2
        for (int j = 0; j < dataPoints; j++){
            myfile2 << arrayM[j] << " " << autocorr[j] << endl;
        }
        // close myfile2
        if (myfile2.is_open()){
            myfile2.close();
        }

        // write in myfile -> tau_e vs temp
        myfile << T[i] << " " << tau_e << endl;

        cout << "T = " << T[i] << ", tau_e = " << tau_e << ", Time: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

    // wrap up
    delete[] arrayM;
    delete[] bootstrap_values;
    delete[] autocorr;
    if (myfile.is_open()){
        myfile.close();
    }
}

void autocorrelation_peak_investigation()
{
    cout << "Running autocorrelation around critical temperature. High computation time! Proceed with caution!" << endl;
    L = 20;
    dim = 2;
    int thermalisationCycles = 15000; // enough time to ensure thermalisation
    int dataPoints = 5000; //"long time", less than before but still long enough
    int spacingCycles = 1;
    double iniT = 2.23; double finT = 2.37; int numT = 8;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,3]:\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,15000);
        cout << "Number of data points [100,20000]\n";
        dataPoints = user_integer_input(100,20000);
        cout << "Initial Temperature [150,400]/100:\n";
        int iniT_int = user_integer_input(150,400);
        iniT = double(iniT_int)/100;
        cout << "Final Temperature [" << iniT_int <<",500]/100:\n";
        int finT_int = user_integer_input(iniT_int,450);
        finT = double(finT_int)/100;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }
    // number of tau_e's to compute and get average
    int numTau;
    cout << "Enter number of tau_e's to compute and average over. Beware of high computation time! [1,1000]:\n";
    numTau = user_integer_input(1,1000);

    //initialise vectors, pointers and time variable
    double *arrayM;
    double *bootstrap_values;
    double *autocorr;
    double *arrayTau;
    int tau_e;
    observable O_tau;
    int t_max = dataPoints;
    clock_t tStartTemp;

    //initialise spins array and neighbours maps
    initialise_system_and_maps();
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Running autocorrelation peak data for L = " << L << endl;

    string folder = ".\\data\\autocorrelation_data\\peak_investigation";
    // open myfile -> average and std
    ofstream myfile;
    string filename = "autocorr_peak_"+to_string(L)+"_L.txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: autocorrelation_peak_investigation(). File not opened with path: " + path + "\nPlease fix path";
    }
    cout << "Writing in file with path: " << path << endl;

    // open myfile2 -> all computed values (to check what is going on)
    ofstream myfile2;
    string filename2 = "autocorr_peak_"+to_string(L)+"_L_allTaus.txt";
    filename_rename_if_exists(filename2,folder);
    string path2 = folder+"\\"+filename2;
    myfile2.open(path2);
    if (!myfile2.is_open()) {
        throw "Func: autocorrelation_peak_investigation(). File not opened with path: " + path2 + "\nPlease fix path";
    }
    cout << "Writing in file with path: " << path2 << endl;

    // begin
    arrayTau = new double[numTau];
    arrayM = new double[dataPoints];
    cout << "Starting computations" << endl;
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        for (int k=0; k<numTau; k++) {

            // compute magnetisation
            initialise_spins_cold();
            compute_magnetisation();
            metropolis_function(T[i],thermalisationCycles);
            arrayM[0] = fabs((double)M/N);
            for (int j=1; j<dataPoints; j++) {
                metropolis_function(T[i],spacingCycles);
                arrayM[j] = fabs((double)M/N);
            }

            // identify tau_e
            autocorr = autocorrelation(arrayM,t_max);
            tau_e = t_max;
            for (int j = 0; j < t_max; j++) {
               if (autocorr[j] <= exp(-1)) {
                   tau_e = j;
                   break;
               }
            }

            // print out, write in myfile2, put in array
            cout << "(" << k+1 << "/" << numTau << ") T = " << T[i] <<  ", tau = " << tau_e << endl;
            myfile2 << T[i] << " " << tau_e << endl;
            arrayTau[k] = double(tau_e);

        }
        //compute average and sigma of arrayTau
        O_tau = compute_average_and_sigma(arrayTau, numTau);

        // write in myfile
        myfile << T[i] << " " << O_tau.value << " " << O_tau.error << endl;
        cout << "\rL = " << L << ", T = " << T[i] << ", tau = " << O_tau.value << " +- " << O_tau.error << ", Time: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

    // wrap up
    delete[] arrayM;
    delete[] arrayTau;
    delete[] bootstrap_values;
    delete[] autocorr;
    if (myfile.is_open()){
        myfile.close();
    }
    if (myfile2.is_open()){
        myfile2.close();
    }
}

/*
ENERGY
Firstly I need create a function to measure the energy of the system. The Hamiltonian H = -J*Sum_<ij>{s_i*s_j} - mu*H*Sum_i{s_i} - {n2n term}. The real question here is if I want to compute the energy manually by calling the function every time I need the E value, or if I want the metropolis_function() to update the energy for me every time (like it does with magnetisation). I don't tink I will need the energy after every MC sweep so I will call to compute it when I need it's value.

Before computing any data I need to manual check that the algorithm iterates over all the edges with no double counting. I will check this with a function, energy_first_check(), and some couts in my compute_energy() function (will comment out after test is passed). TEST PASSED FOR DIM 2 AND 3 :D

Now that the 1st test was passed I will continue without next to nearest neighbours for the time being and conduct a second MANUAL check that the energy calculatation is correct for printed configurations. PASSED

PLOT 1:
Checking how the energy evolves with time at different temperatures and L's is also a sanity check. energy_vs_time_data() will generate my datafiles to be plotted with python

PLOT 2:
plot Energy vs Temperature. I expect to see the the average energy rising from its maximum negative value to 0 as the temperature rises. The average energy is e = E/(nearest_links + next2nearest_links*J') where the denominator is the maximum |energy| the system can achieve. I believe energy_per_link() will be necessary for the plot if I want to plot multiple L's.

*/

void energy_first_check()
{
    cout << "Running for energy of evolving system, check one" << endl;
    dim = 2;
    n2n = 1;
    L = 3;
    int thermalisationCycles = 0;
    int dataPoints = 5;
    int spacingCycles = 10;
    double Temp = 1.5;
    //initialise spins array and neighbours maps
    initialise_system_and_maps();
    //print out nearest and nexttonearest neightbours
    print_mapOfNearest(100);
    print_mapOfNext2Nearest(100);

    initialise_spins_hot();
    //initialise_spins_cold();
    compute_magnetisation();
    print_spins_2D();
    compute_energy();
    cout << "Energy = " << E << ", Mangetisation = " << M << endl;
    for (int i=1; i<dataPoints; i++) {
        metropolis_function(Temp,spacingCycles);
        print_spins_2D();
        compute_energy();
        cout << "Energy = " << E << ", Mangetisation = " << M << endl;
    }
}

void energy_second_check()
{
    cout << "Running for energy of evolving system, check two" << endl;
    dim = 2;
    L = 6;
    int thermalisationCycles = 1000;
    int dataPoints = 10;
    int spacingCycles = 1;
    double Temp = 5;
    //initialise spins array and neighbours maps
    initialise_system_and_maps();
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);
    cout << N_links << endl;
    initialise_spins_hot();
    metropolis_function(Temp,thermalisationCycles);
    //initialise_spins_cold();
    print_spins_2D();
    compute_magnetisation();
    compute_energy();
    cout << "E = " << E << ", e = " << energy_per_link() << ", M = " << M << ", m = " << fabs((double)M/N) << endl;
    for (int i=1; i<dataPoints; i++) {
        metropolis_function(Temp,spacingCycles);
        print_spins_2D();
        compute_energy();
        cout << "E = " << E << ", e = " << energy_per_link() << ", M = " << M << ", m = " << fabs((double)M/N) << endl;
    }
}

void energy_vs_time_data()
{
    cout << "Running for energy vs time (at different temperatures) data" << endl;
    dim = 2;
    L = 80;
    int thermalisationCycles = 0;
    int spacingCycles = 1;
    int dataPoints = 2000;      //total data points
    double iniT = 1.0; double finT = 5.0; int numT = 21;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,4]:\n";
        dim = user_integer_input(2,4);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,10000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    clock_t tStartTemp;
    double e_ave;

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // initiate file
    ofstream myfile;
    string folder = ".\\data\\energy_data\\vs_time";

    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        //open file
        string filename= "energy_data_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(T[i])+".txt";
        filename_rename_if_exists(filename, folder);
        string path = folder+"\\"+filename;
        myfile.open(path);
        if (!myfile.is_open()) {
            throw "Func: energy_vs_temp_data(). File not opened with path: "+path;
        }
        // begin
        initialise_spins_hot();
        compute_magnetisation(); // why not
        metropolis_function(T[i],thermalisationCycles);
        compute_energy();
        myfile << energy_per_link() << endl;
        //arrayE[0] = energy_per_link();
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            compute_energy();
            //arrayE[j] = energy_per_link();
            myfile << energy_per_link() << endl;
        }

        if (myfile.is_open()){
            myfile.close();
        }
    }

   // wrap up
   //delete[] arrayE;
   delete[] spins;
   delete[] T;
}

void energy_vs_temp_data()
{
    cout << "Running for energy at different temperatures data" << endl;
    dim = 2;
    L = 40;
    int thermalisationCycles = 5000;
    int spacingCycles = 50;
    int dataPoints = 2000;      //total data points
    double iniT = 1.0; double finT = 5.0; int numT = 41;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,4]:\n";
        dim = user_integer_input(2,4);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,10000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // pointers and time variable
    double *arrayE;
    double *bootstrap_values;
    clock_t tStartTemp;
    arrayE = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\energy_data\\vs_temp";
    string filename = "energy_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: energy_vs_temp_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // begin
        initialise_spins_auto(T[i]);
        compute_magnetisation();
        metropolis_function(T[i],thermalisationCycles);
        compute_energy();
        arrayE[0] = energy_per_link();
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            compute_energy();
            arrayE[j] = energy_per_link();
        }
        // average and error
        bootstrap_values = bootstrap_error(arrayE, dataPoints, 128, true);

        //write in file and cout things
        myfile << T[i] << " " << bootstrap_values[0] << " " << bootstrap_values[1] << endl;
        cout << "T = " << T[i] << ", E = " << bootstrap_values[0] << " +- " << bootstrap_values[1] << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

   // wrap up
   delete[] arrayE;
   delete[] spins;
   delete[] T;
   if (myfile.is_open()){
       myfile.close();
   }
}

/*
HEAT CAPACITY
Create a function to compute the heat capacity. Initially set it to return c units of the boltzmann constant (hoping that is a good thing to do) to give us a dimensionless heat capacity.

PLOT 1:
Heat capacity vs temperature. Expecting peak close to T_c and peak values increasing with L.

Following that I will compute the heat capacity around the critical temperature. Fit a gaussian on the peak to get the position of the peak T_c(L). Repeat this for multiple L's.

PLOT 2:
T_c(L) vs L and fitted curve to check with Onsager's result.

PLOT 3:
Plot log(T_c(L)-T_c(inf)) vs log(L). Using Onsager's result get values for the constants in the equation with linear regression.
*/

void heat_capacity_data()
{
    cout << "Running for heat capacity at different temperatures data" << endl;
    dim = 2;
    L = 48;
    int thermalisationCycles = 1000;
    int spacingCycles = 50;
    int dataPoints = 2000;      //total data points
    double iniT = 1.0; double finT = 5.0; int numT = 21;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,3]:\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,10000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // pointers and variables
    double *c;
    double *arrayE;
    double *bootstrap_values;
    clock_t tStartTemp;
    arrayE = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\heat";
    string filename = "heat_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: heat_capacity_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // begin
        initialise_spins_auto(T[i]);
        compute_magnetisation(); // why not
        metropolis_function(T[i],thermalisationCycles);
        compute_energy();
        arrayE[0] = E;
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            compute_energy();
            arrayE[j] = E;
        }
        //compute heat capacity
        c = get_heat_capacity(arrayE,dataPoints,T[i]);

        //write in file and cout things
        myfile << T[i] << " " << c[0] << " " << c[1] << endl;
        cout << "T = " << T[i] << ", c = " << c[0] << " +- " << c[1] << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

   // wrap up
   delete[] arrayE;
   delete[] spins;
   delete[] T;
   if (myfile.is_open()){
       myfile.close();
   }
}

void heat_capacity_peak_data()
{
    cout << "Running for heat capacity around peak data" << endl;
    dim = 2;
    L = 20;
    int thermalisationCycles = 1000;
    int spacingCycles = 50;
    int dataPoints = 5000;      //total data points
    double iniT = 2.31; double finT = 2.5; int numT = 20;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,3]:\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,10000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // pointers and variables
    double *c;
    double *arrayE;
    double *bootstrap_values;
    clock_t tStartTemp;
    arrayE = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\heat";
    string filename = "heat_peak_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: heat_capacity_peak_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // begin
        initialise_spins_auto(T[i]);
        compute_magnetisation(); // why not
        metropolis_function(T[i],thermalisationCycles);
        compute_energy();
        arrayE[0] = E;
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            compute_energy();
            arrayE[j] = E;
        }
        //compute heat capacity
        c = get_heat_capacity(arrayE,dataPoints,T[i]);

        //write in file and cout things
        myfile << T[i] << " " << c[0] << " " << c[1] << endl;
        cout << "T = " << T[i] << ", c = " << c[0] << " +- " << c[1] << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

   // wrap up
   delete[] arrayE;
   delete[] spins;
   delete[] T;
   if (myfile.is_open()){
       myfile.close();
   }
}

/*
EXTERNAL FIELD
Investigate what happens with external field: in particular, examine hysteresis effects when H is cycled at different temperatures.

How to do this:
1. Thermalise at H=-H_max. thermalise. record data points. compute average. ouput
2. increment H. thermalise. record data points. compute average. output
3. repeat step 2 until you reach +H_max
4. decrement H. thermalise. record data points. compute average. ouput
5. repeat step 4 until you reach -H_max
stop

This will be done at different Temperatures to investigate how the lattice interacts with the external field below and above T_c.

PLOT 1:
Average Magnetisation vs External Field at different temperatures

PLOT 2:
Avrage Energy vs External Field at different temperatures
*/

void external_field_investigation()
{
    cout << "Running for external field investigation" << endl;
    dim = 2;
    L = 36;
    double H_max = 2.15;
    double H_step = 0.05;
    int thermalisationCycles = 1000;
    int spacingCycles = 50;
    int dataPoints = 100;      //total data points
    double Temp = 1.0;
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,3]:\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,10000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Temperature [1,100]/10\n";
        Temp = double(user_integer_input(1,100))/10;
    }

    // pointers and variables
    bool increase_H = true;
    bool compute = true;
    clock_t tStartTemp;
    double *M_values;
    double *E_values;
    double *arrayE;
    double *arrayM;
    arrayE = new double[dataPoints];
    arrayM = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\H_data";
    string filename = "H_data_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(Temp)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: heat_capacity_peak_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // begin
    H = -H_max;
    initialise_spins_hot();
    compute_magnetisation();

    while (compute) {
        tStartTemp = clock();
        metropolis_function(Temp,thermalisationCycles);
        compute_energy();
        arrayE[0] = E;
        arrayM[0] = M;
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(Temp,spacingCycles);
            compute_energy();
            arrayE[j] = E;
            arrayM[j] = M;
        }
        // average and error
        M_values = bootstrap_error(arrayM, dataPoints, 128, true);
        E_values = bootstrap_error(arrayE, dataPoints, 128, true);
        // write in file
        myfile << H << " " << (M_values[0]) <<  " " << (M_values[1])<< " " <<(E_values[0])  << " " << (E_values[1])<< endl;
        // output
        cout << "H = " << H << ", M = " << M_values[0] <<  " +- " << M_values[1] << ", E = " << E_values[0]  << " +- " << E_values[1] << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
        // de/i-ncrement
        if (increase_H) {
            H += H_step;
            if (H>H_max) {
                increase_H = false;
            }
        } else {
            H -= H_step;
            if (H<-H_max) {
                compute = false;
            }
        }
    }

    if (myfile.is_open()){
        myfile.close();
    }

   // wrap up
   delete[] T;
   delete[] arrayM;
   delete[] arrayE;
   delete[] spins;
}

/*
MAGNETIC SUSCEPTIBILITY
Pretty much the same thing as heat capacity calculations. Will only calculate around peak. Would be a good idea to create a separate analysis function to calculate both chi and C over the same iterations to save computation time. However I already finished computing values for heat capacity therefore will not be necessary.
*/

void magnetic_susceptibility_peak_data()
{
    cout << "Running for magnetic susceptibility around peak data" << endl;
    dim = 2;
    L = 64;
    int thermalisationCycles = 1000;
    int spacingCycles = 50;
    int dataPoints = 10000;      //total data points
    double iniT = 2.20; double finT = 2.40; int numT = 21;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,3]:\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,10000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // pointers and variables
    double *chi;
    double *arrayM;
    double *bootstrap_values;
    clock_t tStartTemp;
    arrayM = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\chi_data";
    string filename = "chi_peak_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: magnetic_susceptibility_peak_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // begin
        initialise_spins_auto(T[i]);
        compute_magnetisation();
        metropolis_function(T[i],thermalisationCycles);
        arrayM[0] = fabs((double)M/N);
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            arrayM[j] = fabs((double)M/N);
        }
        //compute heat capacity
        chi = get_magnetic_susceptibility(arrayM,dataPoints,T[i]);

        //write in file and cout things
        myfile << T[i] << " " << chi[0] << " " << chi[1] << endl;
        cout << "T = " << T[i] << ", chi = " << chi[0] << " +- " << chi[1] << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

   // wrap up
   delete[] arrayM;
   delete[] spins;
   delete[] T;
   if (myfile.is_open()){
       myfile.close();
   }
}

/*
ANIMATION / GIF
Will attempt to create a GIF showing an evolving lattice.
C++ will output the configurations
Python will create the checkerboard pattern of the lattice and will add many successive configurations to create GIFs

Also will create a function that generates a single configuration for a figure
*/

void generate_configurations_for_gif()
{
    cout << "Running to generate configurations for a GIF" << endl;
    dim = 2;
    L = 100;
    int thermalisationCycles = 0;
    int dataPoints = 500;
    int spacingCycles = 1;
    double Temp = 1.0;
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Temperature [1,300]/10\n";
        Temp = double(user_integer_input(1,300))/10;
    }

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\configs\\fig";
    string filename = "configs_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(Temp)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: generate_configurations_for_gif(). File not opened with path: " + path + "\nPlease fix path";
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;

    initialise_spins_hot();
    metropolis_function(Temp,thermalisationCycles);
    // print spins in line
    if (spins[0] == 1) {
        myfile << 1;
    } else {myfile << 0;}
    for (int j=1; j<N; j++) {
        if (spins[j] == 1) {
            myfile << " " << 1;
        } else {myfile << " " << 0;}
    }
    myfile << endl;
    // repeat for dataPoints configurations
    for (int i=1; i<dataPoints; i++) {
        metropolis_function(Temp,spacingCycles);
        // print spins in line
        if (spins[0] == 1) {
            myfile << 1;
        } else {myfile << 0;}
        for (int j=1; j<N; j++) {
            if (spins[j] == 1) {
                myfile << " " << 1;
            } else {myfile << " " << 0;}
        }
        myfile << endl;
    }

    // wrap up
    delete[] spins;
    delete[] T;
    if (myfile.is_open()){
       myfile.close();
    }
}

void generate_configuration_for_figure()
{
    cout << "Running to generate a configuration for a figure" << endl;
    dim = 2;
    L = 150;
    int thermalisationCycles = 3000;
    int dataPoints = 1;
    int spacingCycles = 0;
    double Temp = 1.5;
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Temperature [1,300]/10\n";
        Temp = double(user_integer_input(1,300))/10;
    }

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\configs\\figure";
    string filename = "config_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(Temp)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: generate_configuration_for_figure(). File not opened with path: " + path + "\nPlease fix path";
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;

    initialise_spins_auto(Temp);
    metropolis_function(Temp,thermalisationCycles);
    // print spins in line
    if (spins[0] == 1) {
        myfile << 1;
    } else {myfile << 0;}
    for (int j=1; j<N; j++) {
        if (spins[j] == 1) {
            myfile << " " << 1;
        } else {myfile << " " << 0;}
    }
    myfile << endl;

    // wrap up
    delete[] spins;
    delete[] T;
    if (myfile.is_open()){
       myfile.close();
    }
}

/*
NEXT TO NEAREST NEAREST NEIGHBOURS 2D
Will investigate the effect of introducing next to nearest neigbours interactions. This will be conducted at fixed Lattice size but varying degree of interaction n2n. Since nearest neigbour interaction is set to 1 then the ratio of the interactions is basically the value of n2n.

PLOT 1:
Magnetisation vs Temperature for different n2n

PLOT 2:
Energy per site vs Temperature for different n2n
*/

void next_to_nearest_investigation()
{
    cout << "Running for magnetisation at different temperatures data" << endl;
    dim = 2;
    L = 40;
    n2n = 1.0;
    int thermalisationCycles = 1000;
    int spacingCycles = 50;
    int dataPoints = 2000;      //total data points
    double iniT = 1.0; double finT = 8.0; int numT = 71;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,4]:\n";
        dim = user_integer_input(2,4);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Next to nearest interaction [0,200]/100:\n";
        n2n = user_integer_input(0,200);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,41]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // initialise pointers and time variable
    double *arrayE;
    double *arrayM;
    double *bootstrap_values_M;
    double *bootstrap_values_E;
    clock_t tStartTemp;
    arrayM = new double[dataPoints];
    arrayE = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\n2n_data";
    string filename = "n2n_data_"+to_string(dim)+"D_"+to_string(L)+"_"+to_string(n2n)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: next_to_nearest_investigation(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    // compute data
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // begin
        initialise_spins_auto(T[i]/(1.0+n2n)); // qualitative. basically increasing n2n increases T_c therefore it would be more efficient to cold start at higher temperatures
        compute_magnetisation();
        metropolis_function(T[i],thermalisationCycles);
        arrayM[0] = fabs((double)M/N);
        arrayE[0] = energy_per_site();
        for (int j=1; j<dataPoints; j++) {
            metropolis_function(T[i],spacingCycles);
            arrayM[j] = fabs((double)M/N);
            arrayE[j] = energy_per_site();
        }
        // average and error
        bootstrap_values_M = bootstrap_error(arrayM, dataPoints, 128, true);
        bootstrap_values_E = bootstrap_error(arrayE, dataPoints, 128, true);
        // write data in myfile
        myfile << T[i] << " " << bootstrap_values_M[0] << " " << bootstrap_values_M[1] << " " << bootstrap_values_E[0] << " " << bootstrap_values_E[1] << endl;
        cout << "T = " << T[i] << ", m = " << bootstrap_values_M[0] << " +- " <<  bootstrap_values_M[1] << ", e = " << bootstrap_values_E[0] << " +- " <<  bootstrap_values_E[1] << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

   // wrap up
   delete[] spins;
   delete[] T;
   delete[] arrayM;
   delete[] bootstrap_values_M;
   delete[] arrayE;
   delete[] bootstrap_values_E;
   if (myfile.is_open()){
       myfile.close();
   }
}

/*
To overcome critical slowing down we will implement the Wolff algorithm. The major difference with this algorithm is that it constructs clusters of sites to collectively flip their spin.

The main difference is noticeable at T_c where large clusters of like spin tend to form. The Metropolis reforms these clusters by slowly adjusting their boundaries, leading to high autocorrelations and therefore the critical slowing down effect. It is virtually impossible for a cluster of opposite spin to be formed inside these clusters. The Wolff has the ability to reform these clusters fast and thus does not suffer from critical slowing down. At high T the Wolff clusters are close to 1 site and the systems evolution resembles closely that of Metropolis

Firsly, we will examine the typical cluster size of the algorithm over a range of temperatures and at different lattice sizes.

Secondly I will examine the autocorrelation of the Wolff in a range around the critical temperature to qualitatevly compare with Metropolis. I expect tau_e to be of order unity throughout

PLOT 1:
Cluster size vs Temperature

PLOT 2:
Tau_e vs Temperature around the critical point

*/

void wolff_cluster_size_vs_temp_data()
{
    cout << "Running for cluster size at different temperatures data" << endl;
    dim = 2;
    L = 30;
    int thermalisationCycles = 100;
    int spacingCycles = 1;
    int dataPoints = 1000;      //total data points
    double iniT = 2.425; double finT = 2.425; int numT = 1;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,4]:\n";
        dim = user_integer_input(2,4);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Next to nearest interaction [0,200]/100:\n";
        n2n = user_integer_input(0,200);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points [100,10000]:\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature [1,99]/10:\n";
        int iniT_int = user_integer_input(1,99);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature [" << iniT_int <<",100]/10:\n";
        int finT_int = user_integer_input(iniT_int,100);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,41]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // initialise pointers and time variable
    double *n; // cluster size average (no error)
    double *array_n;
    array_n = new double[dataPoints];
    clock_t tStartTemp;
    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open file
    ofstream myfile;
    string folder = ".\\data\\wolff_data\\cluster_size";
    string filename = "clustersize_data_"+to_string(dim)+"D_"+to_string(L)+".txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: wolff_cluster_size_vs_temp_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;

    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        // begin
        initialise_spins_auto(T[i]);
        compute_magnetisation();
        n = wolff_function(T[i],thermalisationCycles);
        for (int j=0; j<dataPoints; j++) {
            n = wolff_function(T[i],spacingCycles);
            array_n[j] = (*n)/N;
        }
        // average and error
        observable n_values = compute_average_and_sigma(array_n, dataPoints);
        // write data in myfile
        myfile << T[i] << " " << n_values.value << " " << n_values.error << endl;
        cout << "T = " << T[i] << ", <n> = " << n_values.value << " +- " <<  n_values.error << ", Time taken: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

    // wrap up
    delete[] spins;
    delete[] T;
    delete[] array_n;
    if (myfile.is_open()){
        myfile.close();
    }
}

void wolff_autocorrelation_investigation()
{
    cout << "Running autocorrelation around critical temperature. High computation time! Proceed with caution!" << endl;
    L = 30;
    dim = 2;
    int thermalisationCycles = 1000; // enough time to ensure thermalisation
    int dataPoints = 5000; //"long time"
    int spacingCycles = 1;
    double iniT = 2.6; double finT = 3.5; int numT = 9;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 1 for YES, 0 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 0) {
        cout << "Dimensions [2,3]:\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size [8,181]:\n";
        L = user_integer_input(8,181);
        cout << "Thermalisation cycles [0,10000]:\n";
        thermalisationCycles = user_integer_input(0,15000);
        cout << "Number of data points [100,20000]\n";
        dataPoints = user_integer_input(100,20000);
        cout << "Initial Temperature [150,400]/100:\n";
        int iniT_int = user_integer_input(150,400);
        iniT = double(iniT_int)/100;
        cout << "Final Temperature [" << iniT_int <<",500]/100:\n";
        int finT_int = user_integer_input(iniT_int,450);
        finT = double(finT_int)/100;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points [2,100]:\n";
            numT = user_integer_input(2,100);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }
    // number of tau_e's to compute and get average
    int numTau;
    cout << "Enter number of tau_e's to compute and average over. Beware of high computation time! [1,1000]:\n";
    numTau = user_integer_input(1,1000);

    //initialise vectors, pointers and time variable
    double *wolff;
    double *arrayM;
    double *bootstrap_values;
    double *autocorr;
    double *arrayTau;
    int tau_e;
    observable O_tau;
    clock_t tStartTemp;
    arrayTau = new double[numTau];
    arrayM = new double[dataPoints];

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    string folder = ".\\data\\wolff_data\\autocorrelation";
    // open myfile -> average and std
    ofstream myfile;
    string filename = "wolff_autocorr_"+to_string(L)+"_L.txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: wolff_autocorrelation_investigation(). File not opened with path: " + path + "\nPlease fix path";
    }
    cout << "Writing in file with path: " << path << endl;

    // begin
    cout << "Starting computations" << endl;
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        for (int k=0; k<numTau; k++) {
            // compute magnetisation
            initialise_spins_auto(T[i]);
            compute_magnetisation();
            wolff = wolff_function(T[i],thermalisationCycles);
            arrayM[0] = fabs((double)M/N);
            for (int j=1; j<dataPoints; j++) {
                wolff = wolff_function(T[i],spacingCycles);
                arrayM[j] = fabs((double)M/N);
            }

            // identify tau_e
            autocorr = autocorrelation(arrayM,dataPoints);
            tau_e = dataPoints;
            for (int j = 0; j < dataPoints; j++) {
               if (autocorr[j] <= exp(-1)) {
                   tau_e = j;
                   break;
               }
            }

            // print out, put in array
            cout << "(" << k+1 << "/" << numTau << ") T = " << T[i] <<  ", tau = " << tau_e << endl;
            arrayTau[k] = double(tau_e);
        }
        //compute average and sigma of arrayTau
        O_tau = compute_average_and_sigma(arrayTau, numTau);

        // write in myfile
        myfile << T[i] << " " << O_tau.value << " " << O_tau.error << endl;
        cout << "\rL = " << L << ", T = " << T[i] << ", tau = " << O_tau.value << " +- " << O_tau.error << ", Time: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

    // wrap up
    delete[] arrayM;
    delete[] arrayTau;
    delete[] bootstrap_values;
    delete[] autocorr;
    if (myfile.is_open()){
        myfile.close();
    }
}
