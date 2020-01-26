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
TESTING THE PROGRAM'S CORE FUNCTIONS
To start with analysing the evolution of the average magnetisation of the system. I will do this with different Lattice sizes and a range of temperatures. Following that I will test how the average magnetisation of an evolved system changes with temperature. I will do this for different lattice sizes to show finite size effects.

PLOT 1:
For a fixed lattice size plot the evolution of the system at different temperatures, both below and above Tc.

PLOT 2:
I will choose a specific lattice size and temperature (below Tc) to demonstrate the randomness of the simulation and difference between hot and cold start

PLOT 3:
For a fixed temperature (below Tc) plot the evolution of the system at different lattice sizes. Use average magnetisation for clarity in plots.

PLOT 4:
<|magnetisaiton|> vs temperature for different lattice sizes. Shows finite size effects and how greater lattice size are closer to the real solution.
*/

void magnetisation_vs_time_data()
{
    cout << "Running for magnetisation of evolving system" << endl;
    dim = 2;
    L = 40;
    int thermalisationCycles = 0;
    int dataPoints = 5000;
    int spacingCycles = 1;
    double Temp = 1.5;
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, 0, Temp);
    cout << "Proceed with default parameters? Enter 0 for YES, 1 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 1) {
        cout << "Dimensions (2-3):\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size (8-256):\n";
        L = user_integer_input(8,256);
        cout << "Thermalisation cycles (0-5000):\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points (100-10000):\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Temperature (enter 15 for 1.5 Kelvin) (1-50):\n";
        Temp = double(user_integer_input(1,50))/10;
    }

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

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
    cout << "Running for magnetisation of evolving systems" << endl;
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
    cout << "Running for magnetisation at different temperatures data" << endl;
    dim = 2;
    L = 64;
    int thermalisationCycles = 10000;
    int spacingCycles = 50;
    int dataPoints = 2000;      //total data points
    double iniT = 1.0; double finT = 5.0; int numT = 41;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 0 for YES, 1 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 1) {
        cout << "Dimensions (2-3):\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size (8-256):\n";
        L = user_integer_input(8,256);
        cout << "Thermalisation cycles (0-5000):\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points (100-10000):\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature (1-49) i.e enter 15 for 1.5 Kelvin:\n";
        int iniT_int = user_integer_input(1,49);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature (" << iniT_int <<"-50) i.e enter 15 for 1.5 Kelvin:\n";
        int finT_int = user_integer_input(iniT_int,50);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points (2-41).\n";
            numT = user_integer_input(2,41);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    // initialise vectors, pointers and time variable
    vector<double> magn;
    vector<double> err_magn;
    double *arrayM;
    double *bootstrap_values;
    clock_t tStartTemp;

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

    // write data in myfile
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

/*
AUTOCORRELATION
1. Firstly, autocorrelation_investigation() will investigate the autcorrelation function and it's behaviour over a wide range of temperatures. A single value will be recorded at each temperature where the autocorr drops by 1/e, tau_e.
2. Lastly, we want to examine the value of tau_e around the critical temperature. From the previous analysis it is evident that there is a large variance over the value even with the same parameters.

N.B.:
It was observed that the total number of data points you collect for your autocorrelation analysis noticably affects the values of tau_e around the critical point. More data points -> higher values of tau_e. Without wasting (computational) time investigating this effect I fixed the number of data points collected to 1000, a small number to reduce computational time and big enough to yield  data for bigger L's. The reduced computational time can be used for more tau_e calculations with the same parameters to get a more precise average.


*/

void autocorrelation_investigation()
{
    cout << "Running autocorrelation investigation" << endl;
    L = 32;
    dim = 2;
    int thermalisationCycles = 3000;
    int dataPoints = 1000;
    int spacingCycles = 1;
    double iniT = 1.0; double finT = 5.0; int numT = 11;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 0 for YES, 1 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 1) {
        cout << "Dimensions (2-3):\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size (8-256):\n";
        L = user_integer_input(8,256);
        cout << "Thermalisation cycles (0-5000):\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points (100-10000):\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature i.e enter 15 for 1.5 Kelvin (1-49):\n";
        int iniT_int = user_integer_input(1,49);
        iniT = double(iniT_int)/10;
        cout << "Final Temperature i.e enter 15 for 1.5 Kelvin (" << iniT_int <<"-50):\n";
        int finT_int = user_integer_input(iniT_int,50);
        finT = double(finT_int)/10;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points (2-41).\n";
            numT = user_integer_input(2,41);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }

    //initialise vectors, pointers and time variable
    double *arrayM;
    double *bootstrap_values;
    clock_t tStartTemp;
    double *autocorr;
    int tau_e;

    //initialise spins array and neighbours maps
    initialise_system_and_maps();

    // open myfile
    ofstream myfile;
    ofstream myfile2;
    string folder = ".\\data\\autocorrelation_data\\investigation";
    string filename = "autocorr_times_"+to_string(L)+"_L.txt";
    string filename2;
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    string path2;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: autocorrelation_investigation(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    for (int i = 0; i<numT; i++) {
        // open myfile2
        filename2 =  "autocorr_data_"+to_string(L)+"_L_"+to_string(T[i])+"_T.txt";
        filename_rename_if_exists(filename2, folder);
        path2 = folder+"\\"+filename2;
        myfile2.open(path2);
        if (!myfile2.is_open()) {
            throw "Func: autocorrelation_investigation() (in temperatures for-loop). File not opened with path: " + path2;
        }
        tStartTemp = clock();
        arrayM = new double[dataPoints];
        // begin
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
        tau_e = 0;
        for (int j = 0; j < dataPoints; j++){
           if (autocorr[j] <= exp(-1)){
               tau_e = j;
               break;
           }
        }

        // write in myfile
        myfile << T[i] << " " << tau_e << endl;

        // write in myfile2
        for(int j=0; j<dataPoints; j++) {
           myfile2 << arrayM[j] << " " << autocorr[j] << endl;
        }
        if (myfile2.is_open()){
            myfile2.close();
        }

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

void autocorrelation_peaks_data()
{
    cout << "Running autocorrelation around critical tepmerature" << endl;
    L = 28;
    dim = 2;
    int thermalisationCycles = 1000;
    int dataPoints = 1000;
    int spacingCycles = 1;
    double iniT = 2.21; double finT = 2.39; int numT = 10;
    T = linspace(iniT, finT, numT);
    print_all_parameters(thermalisationCycles, dataPoints, spacingCycles, numT, 0);
    cout << "Proceed with default parameters? Enter 0 for YES, 1 for NO\n";
    int user_input = user_integer_input(0,1);
    if (user_input == 1) {
        cout << "Dimensions (2-3):\n";
        dim = user_integer_input(2,3);
        cout << "Lattice size (8-256):\n";
        L = user_integer_input(8,256);
        cout << "Thermalisation cycles (0-5000):\n";
        thermalisationCycles = user_integer_input(0,5000);
        cout << "Number of data points (100-10000):\n";
        dataPoints = user_integer_input(100,10000);
        cout << "Initial Temperature i.e enter 150 for 1.50 Kelvin (10-490):\n";
        int iniT_int = user_integer_input(10,490);
        iniT = double(iniT_int)/100;
        cout << "Final Temperature i.e enter 150 for 1.50 Kelvin (" << iniT_int <<"-500):\n";
        int finT_int = user_integer_input(iniT_int,500);
        finT = double(finT_int)/100;
        if (iniT_int != finT_int) {
            cout << "Number of Temperature points (2-41).\n";
            numT = user_integer_input(2,41);
        }
        else {
            numT = 1;
        }
        T = linspace(iniT, finT, numT);
    }
    // number of taus to compute and get average
    int numTau; // good number is 1000
    cout << "Enter number of tau_e's to compute:\n";
    numTau = user_integer_input(10,1000);

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

    // open myfile
    ofstream myfile;
    string folder = ".\\data\\autocorrelation_data\\peaks";
    string filename = "autocorr_peak_"+to_string(L)+"_L.txt";
    filename_rename_if_exists(filename, folder);
    string path = folder+"\\"+filename;
    myfile.open(path);
    if (!myfile.is_open()) {
        throw "Func: autocorrelation_peaks_data(). File not opened with path: "+path;
    }
    cout << "Writing in file with path: " << path << endl;
    cout << "Starting computations" << endl;
    for (int i = 0; i<numT; i++) {
        tStartTemp = clock();
        arrayTau = new double[numTau];
        for (int k=0; k<numTau; k++) {
            cout << "\r" << k*100/numTau << " %  ";
            arrayM = new double[dataPoints];
            // begin
            initialise_spins_cold();
            compute_magnetisation();
            metropolis_function(T[i],thermalisationCycles);
            arrayM[0] = M;
            for (int j=1; j<dataPoints; j++) {
                metropolis_function(T[i],spacingCycles);
                arrayM[j] = M;
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
            arrayTau[k] = double(tau_e);
        }
        //compute average and sigma of arrayTau
        O_tau = compute_average_and_sigma(arrayTau, numTau);

        // write in myfile
        myfile << T[i] << " " << O_tau.value << " " << O_tau.error << endl;
        cout << "\rL = " << L << ", T = " << T[i] << ", tau = " << O_tau.value << " +- " << O_tau.error << ", Time: " << (double)(clock()-tStartTemp)/CLOCKS_PER_SEC << " seconds" << endl;
    }

    // wrap up
    delete[] arrayTau;
    delete[] arrayM;
    delete[] bootstrap_values;
    delete[] autocorr;
    if (myfile.is_open()){
        myfile.close();
    }
}

int initial_menu()
{
    cout << "Ising Model, made by Savvas Shiakas (ss2477)" << endl;
    cout << "C++ was used to produce data files which are later plotted using python\n";
    cout << "Please select an analysis by entering an integer:" << endl;
    cout << "1 for Magnetisation vs Time.\n";
    cout << "2 for Magnetisation vs Time (Bulk data, different lattice sizes and temperatures).\n";
    cout << "3 for Magnetisation vs Temperature.\n";
    cout << "4 for Autocorrelation initial Investigation.\n";
    cout << "5 for Autocorrelation vs Temperature around the critical temperature - high computation time.\n";
    cout << "0 to Exit the program.\n";
    int choice = user_integer_input(0,5);
    return choice;
}

int main(int argc, char** argv)
{
    clock_t tStart = clock();
    srand(seed);

    // TODO: enclose in while loop for code to rerun until user selects 0
    // user friendly function to run the code
    bool running_program = true;
    while (running_program) {
        int user_choice = initial_menu();
        try {
            switch(user_choice) {
                case 0: running_program = false;
                break;
                case 1: magnetisation_vs_time_data();
                break;
                case 2: magnetisation_vs_time_data_bulk();
                break;
                case 3: magnetisation_vs_temp_data();
                break;
                case 4: autocorrelation_investigation();
                break;
                case 5: autocorrelation_peaks_data();
                break;
            }
        }
        catch(string er) {
            cerr << "ERROR. EXITING PROGRAM." << endl;
            cout << er << endl;
            return 1;
        }
    }

    cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    return 0;
}
