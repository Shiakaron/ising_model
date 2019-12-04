//library
#include "header.h"

//GLOBAL VARIABLES
int L = 32;
int N = L*L;                //number of sites
int **F;                    //spin array pointer (2dim)
int *F_2;                   //spin array pointer (Ndim)
unsigned int seed =  (unsigned)time(0); // (unsigned)time(0);

int nCycles = 10000;        //number of cycles of Monte Carlo primary
int mCycles = 50;           //number of cycles of Monte Carlo secondary
int dataPoints = 1000;      //total data points

int n_b = 8;                // n_b is number of bins, used in jackknife error analysis

double iniT = 2.0;          //initial temperature
double finT = 7.0;          //final temperature
int numT = 51;              //data points (for linear spacing)
double *T;                  //array of temperatures

double tempM;               //used for less time consuming magnetization calculation

double Next2Nearest = 0;
double B = 0;
int dim = 2;
int V = pow(L,dim);

void MagnetisationEnergyData() {
    /*
     Function that computes data to plot Magnetization and Energy vs Temperature
     */
    // initialize temperature array
    T = linspace(iniT, finT, numT);

    ofstream myfile; // initialiZe file to input data points
    myfile.open("data.txt"); // TXT file with columns: Temperature, Magnetization, Energy, SigmaJackM, SigmaJackE


    double aveM; // initialize the variable names for the average Energy and Magnetization at each Temperature
    double aveE;

    double sigmaJackM;
    double sigmaJackE;

    double arrayM[dataPoints]; // initialize reusable arrays to insert data points, these arrays are used to compute an average and a sigma at each Temperature.
    double arrayE[dataPoints]; // the +1 is for the 1st data point after nCycles and the rest is for the secondary phase with mCycles in between them

    // begin for-loop to compute values for each Temperature.
    for (int i =0; i<numT; i++) {
        // Loading
        float pct = 100*i/numT;
        cout << "\r" << "Data Loading : " << setprecision(2) << pct << "%";

        myfile << T[i]<< ","; // 1st column in data.txt to be the Temperature

        initiStateHot(); // initialize the state
        flipFunction(T[i], nCycles); // perform nCycles of the Monte Carlo sweep (Metropolis)


        // insert the 1st data points after nCycles
        arrayM[0] = fabs(computeMagnetisation());
        arrayE[0] = computeEnergy();

        // for loop for secondary phase to get an average for Energy and Magnetization

        for (int j = 1; j < dataPoints; j++){
            flipFunction(T[i],mCycles);

            arrayM[j]= fabs(computeMagnetisation());
            arrayE[j]=computeEnergy();
        }


        //compute averages
        aveM = averageArray(arrayM, dataPoints);
        aveE = averageArray(arrayE, dataPoints);

        //compute jackknife errors (better error estimate)
        sigmaJackM = jackknife_new(arrayM, n_b);
        sigmaJackE = jackknife_new(arrayE, n_b);

        //print out the averages and sigmas on the TXT file
        myfile << aveM << ",";
        myfile << aveE << ",";
        myfile << sigmaJackM << ",";
        myfile << sigmaJackE << "\n";

    }

    //end of sequence, close file and print out that the computations are finished.
    cout << "\r" << "Data Loading : 100%" << endl;
    if(myfile.is_open()){
        myfile.close();
        cout << "finished data" << endl;
    }

}

void ThermalisationInvestigation() {
    /*
     Thermalization investigation. The function prints out in a file the magnetisation of the sample on each step of
     the Metropolis algorithm.
     */
    ofstream myfile5;
    myfile5.open("Thermalisation.txt"); // Columns: sweep, <Magnetisation>

    //Temperature
    double T = 2.2;
    //initializations

    int multi = 1;
    int xCycles = multi*nCycles;
    double M = 0;
    initiStateHot();

    for (int i = 0; i<xCycles; i++) {
        // Loading
        //float pct = 100*i/xCycles;
        //cout << "\r" << "Thermalisation Loading : " << setprecision(2) << pct << "%";

        M = fabs(computeMagnetisation());

        myfile5 << i << "," << M << endl;
        flipFunction(T, 1);
    }

    if(myfile5.is_open()) {
        myfile5.close();
        cout << "finished Thermalisation" << endl;
    }

    cout << "Final magnetisation: " << M << endl;
}

void autocorrelationInvestigation(){
    /*
    The autocorrelation investigation function as a first stage acts exactly like he thermalisation investigation but stores
    the magnetisation as a function of cycles in an array instead of printing it off on a file. Then the autocorrelation
    function takes over to compute the function's value which is later printed out along with magnetisation.
    */
    //initiate TXT file
    ofstream myfile6;
    myfile6.open("Autocorrelation.txt");
    //initialise temperature and magnetisation array and autocorrelation pointer
    double T = 2.2;
    double M[nCycles];
    double *autocorr;

    //initialise state
    initiStateHot();


    M[0] = fabs(computeMagnetisation());

    //
    for (int i = 1; i<nCycles; i++) {
        // Loading
        float pct = 100*i/nCycles;
        cout << "\r" << "Autocorrelation Loading : " << setprecision(2) << pct << "%";

        flipFunction(T, 1);
        M[i] = fabs(computeMagnetisation());
    }
    cout << "\r" << "Autocorrelation Loading : 100%" << endl;
    cout << "just one moment please " << endl;

    autocorr = autocorrelation(M,nCycles);

    for (int i = 0; i<nCycles; i++) {
        myfile6 << i << "," << M[i] << "," << autocorr[i] << endl;
    }

    //end of sequence, close file and print out that the computations are finished.
    if(myfile6.is_open()) {
        myfile6.close();
        cout << "finished autocorrelation" << endl;
    }
}

/************ Same program using more less time consuming Magnetization calculation ***********/

void MagnetisationEnergyDataNew() {
    /*
     Function that computes data to plot Magnetization and Energy vs Temperature
     */
    // initialize temperature array
    T = linspace(iniT, finT, numT);

    ofstream myfile; // initialiZe file to input data points
    myfile.open("dataN.txt"); // TXT file with columns: Temperature, Magnetization, Energy, SigmaJackM, SigmaJackE

    double aveE; // initialize the variable names for the average Energy and Magnetization at each Temperature
    double aveMN;

    double sigmaJackE;
    double sigmaJackMN;

    double arrayE[dataPoints]; // the +1 is for the 1st data point after nCycles and the rest is for the secondary phase with mCycles in between them
    double arrayMNEW[dataPoints]; // initialize reusable arrays to insert data points, these arrays are used to compute an average and a sigma at each Temperature.

    // begin for-loop to compute values for each Temperature.
    for (int i =0; i<numT; i++) {
        // Loading
        float pct = 100*i/numT;
        cout << "\r" << "Data Loading : " << setprecision(2) << pct << "%";

        myfile << T[i]<< ","; // 1st column in data.txt to be the Temperature

        initiStateHot(); // initialize the state
        flipFunction(T[i], nCycles); // perform nCycles of the Monte Carlo sweep (Metropolis)
        tempM = computeMagnetisation();

        // insert the 1st data points after nCycles
        arrayMNEW[0] = abs(tempM);
        arrayE[0] = computeEnergy();

        // for loop for secondary phase to get an average for Energy and Magnetization
        for (int j = 1; j < dataPoints; j++){
            flipFunction(T[i],mCycles);
            arrayMNEW[j] = fabs(tempM);
            arrayE[j]=computeEnergy();
        }

        //compute averages
        aveE = averageArray(arrayE, dataPoints);
        aveMN = averageArray(arrayMNEW, dataPoints);

        //compute jackknife errors (better error estimate)
        sigmaJackE = jackknife_new(arrayE, n_b);
        sigmaJackMN = jackknife_new(arrayMNEW, n_b);

        //print out the averages and sigmas on the TXT file
        myfile << aveMN << ",";
        myfile << aveE << ",";
        myfile << sigmaJackMN << ",";
        myfile << sigmaJackE << "\n";
    }

    //end of sequence, close file and print out that the computations are finished.
    cout << "\r" << "Data Loading : 100%" << endl;
    if(myfile.is_open()){
        myfile.close();
        cout << "finished data" << endl;
    }

}

void ThermalisationInvestigationNew() {
    /*
     Thermalization investigation. The function prints out in a file the magnetisation of the sample on each step of
     the Metropolis algorithm.
     */
    ofstream myfile5;
    myfile5.open("ThermalisationN.txt"); // Columns: sweep, <Magnetisation>

    //Temperature
    double T = 2.3;

    //initializations
    int multi = 1;
    int xCycles = multi*nCycles;
    initiStateHot();
    tempM = (computeMagnetisation());

    for (int i = 0; i<xCycles; i++) {
        // Loading
        float pct = 100*i/xCycles;
        cout << "\r" << "ThermalisationNew Loading : " << setprecision(2) << pct << "%";

        //M = computeMagnetisation();
        myfile5 << i << "," << fabs(tempM) << endl;
        flipFunction(T, 1);
    }
    cout << "\r" << "ThermalisationNew Loading : 100%" << endl;
    if(myfile5.is_open()) {
        myfile5.close();
        cout << "finished ThermalisationNew" << endl;
    }
}

void autocorrelationInvestigationNew(){
    /*
    The autocorrelation investigation function as a first stage acts exactly like he thermalisation investigation but stores
    the magnetisation as a function of cycles in an array instead of printing it off on a file. Then the autocorrelation
    function takes over to compute the function's value which is later printed out along with magnetisation.
    */
    //initiate TXT file
    ofstream myfile6;
    myfile6.open("AutocorrelationN.txt");
    //initialize temperature and magnetization array and autocorrelation pointer
    double T = 2.2;
    double M[nCycles];
    double *autocorr;

    //initialize state
    initiStateHot();

    tempM = (computeMagnetisation());
    M[0] = abs(tempM);

    //
    for (int i = 1; i<nCycles; i++) {
        // Loading
        float pct = 100*i/nCycles;
        cout << "\r" << "Autocorrelation Loading : " << setprecision(2) << pct << "%";

        flipFunction(T, 1);

        M[i] = fabs(tempM);
    }
    cout << "\r" << "Autocorrelation Loading : 100%" << endl;
    cout << "just one moment please " << endl;

    autocorr = autocorrelation(M,nCycles);

    for (int i = 0; i<nCycles; i++) {
        myfile6 << i << "," << M[i] << "," << autocorr[i] << endl;
    }

    //end of sequence, close file and print out that the computations are finished.
    if(myfile6.is_open()) {
        myfile6.close();
        cout << "finished autocorrelation" << endl;
    }
}

/*********** Another spin off (pun intended) of the MagnetisationEnergyData() function to simulate spin configurations *************/

void MagnetisationEnergyData2() {
    /*
     Function that computes data to plot Magnetisation and Energy vs Temperature
     Added spins.txt to print out spins and then create the visual simulation in matlab
     */
    // initialise temperature array
    T = linspace(iniT, finT, numT);

    ofstream myfile; // initialise file to input data points
    myfile.open("data.txt"); // TXT file with columns: Temperature, Magnetisation, Energy, SigmaJackM, SigmaJackE

    ofstream myfile2;
    myfile2.open("spins.txt");

    double aveM; // initialise the variable names for the average Energy and Magnetisation at each Temperature
    double aveE;

    double sigmaJackM;
    double sigmaJackE;

    double arrayM[dataPoints]; // initialise reusable arrays to insert data points, these arrays are used to compute an average and a sigma at each Temperature.
    double arrayE[dataPoints]; // the +1 is for the 1st data point after nCycles and the rest is for the secondary phase with mCycles in between them

    // begin for-loop to compute values for each Temperature.
    for (int i =0; i<numT; i++) {
        // Loading
        float pct = 100*i/numT;
        cout << "\r" << "Data Loading : " << setprecision(2) << pct << "%";

        myfile << T[i]<< ","; // 1st column in data.txt to be the Temperature

        initiStateHot(); // initialise the state
        flipFunction(T[i], nCycles); // perform nCycles of the Monte Carlo sweep (Metropolis)

        // insert the 1st data points after nCycles
        arrayM[0] = computeMagnetisation();
        arrayE[0] = computeEnergy();

        bool printToSpinsFile = (T[i]==1.8); // temperature to print spins

        if (printToSpinsFile) {printArraytoFile(myfile2);}

        // for loop for secondary phase to get an average for Energy and Magnetisation
        for (int j = 1; j < dataPoints; j++){
            flipFunction(T[i],mCycles);
            arrayM[j]=computeMagnetisation();
            arrayE[j]=computeEnergy();
            if (printToSpinsFile) {printArraytoFile(myfile2);}
        }

        //compute averages
        aveM = averageArray(arrayM, dataPoints);
        aveE = averageArray(arrayE, dataPoints);

        //compute jackknife errors (better error estimate)
        sigmaJackM = jackknife_new(arrayM, n_b);
        sigmaJackE = jackknife_new(arrayE, n_b);

        //print out the averages and sigmas on the TXT file
        myfile << aveM << ",";
        myfile << aveE << ",";
        myfile << sigmaJackM << ",";
        myfile << sigmaJackE << "\n";

    }

    //end of sequence, close file and print out that the computations are finished.
    cout << "\r" << "Data Loading : 100%" << endl;
    if(myfile2.is_open()){
        myfile2.close();
    }
    if(myfile.is_open()){
        myfile.close();
        cout << "finished data" << endl;
    }

}

void simulationOfSpins(double T, int sCycles, int nConfigurations) {
    /*
    This function will output on a text file nConfigurations of spins
    with sCycles of spacing between them, which can be later be made
    into a GIF using matlab. This will work for a single temperature at
    a time.
    */
    ofstream myfile2;
    myfile2.open("spins.txt");

    initiStateHot(); // initialise the state
    printArraytoFile(myfile2);

    // for loop for secondary phase to get an average for Energy and Magnetisation
    for (int j = 1; j < nConfigurations; j++){
        flipFunction(T,sCycles);
        printArraytoFile(myfile2);
    }

    //end of sequence, close file and print out that the computations are finished.
    if(myfile2.is_open()){
        myfile2.close();
        cout << "spins finished " << endl;
    }

}

void simulationOfSpins2(double T, int sCycles, int nConfigurations, bool thermalised) {
    /*
    Final version for the simulation of spins. The user needs to input the temperature,
    the cycles in between configurations, the number of configurations he wants to be
    printed out and if the state is already thermalised at the 1st configuration. The
    printed file also holds this information on its name. The style of printing out is a
    bit different as well, instead of printing out the matrix as in the previous version we
    have the first 2 lines describing the parameters of the configurations and then 3
    columns giving the x and y coordinates of the spin, up/down.
    */
    ofstream myfile;

    if (thermalised == true) {
        myfile.open("spins_thermalised_"+to_string(L)+"_"+to_string(T)+"_"+to_string(sCycles)+"_"+to_string(nConfigurations)+".txt");
        initiStateHot(); // initialise the state
        flipFunction(T,nCycles); // thermalise
        //myfile << "Lattice size ,Temperature , Cycles in between,  Configurations computed, thermalised(1/0)"<< endl;
        //myfile << L << "," << T << "," << sCycles << "," << nConfigurations << "," << thermalised << endl;
        //myfile << "x,y,up/down" << endl;
        printArraytoFile2(myfile); // print 1st configuration

        for (int j = 1; j < nConfigurations; j++){
            flipFunction(T,sCycles);
            printArraytoFile2(myfile);
        }

        //end of sequence, close file and print out that the computations are finished.
        if(myfile.is_open()){
            myfile.close();
            cout << "spins_thermalised_"+to_string(L)+"_"+to_string(T)+"_"+to_string(sCycles)+"_"+to_string(nConfigurations)+" finished " << endl;
        }

    }
    else {
        myfile.open("spins_unthermalised_"+to_string(L)+"_"+to_string(T)+"_"+to_string(sCycles)+"_"+to_string(nConfigurations)+".txt");
        initiStateHot(); // initialise the state
        //no thermalisation taking place
        myfile << "Lattice size ,Temperature , Cycles in between,  Configurations computed, thermalised(1/0)"<< endl;
        myfile << L << "," << T << "," << sCycles << "," << nConfigurations << "," << thermalised << endl;
        myfile << "x,y,up/down" << endl;
        printArraytoFile2(myfile); // print 1st configuration

        for (int j = 1; j < nConfigurations; j++){
            flipFunction(T,sCycles);
            printArraytoFile2(myfile);
        }

        //end of sequence, close file and print out that the computations are finished.
        if(myfile.is_open()){
            myfile.close();
            cout << "spins_unthermalised_"+to_string(L)+"_"+to_string(T)+"_"+to_string(sCycles)+"_"+to_string(nConfigurations)+" finished " << endl;
        }
    }
}

/***********  computation of the susceptibility of the lattice ***********/

void susceptibilityData() {
    /*
    Investigating the susceptibility, chi, of the sample over a range of temperatures.
    This will give out data for various lattice sizes L around the critical temperature.
    Since we only really care about the region around T_crit we will set iniT = 2 and finT = 3
    and compute enough points so that fitting a curve will not be necessary.
    We will extract only the maximum value of chi for each L.
    */

    //initialise the files we will be writing in
    ofstream myfile;
    myfile.open("Tpseudo.txt"); // Columns: L, T_pseudo, T_error

    ofstream myfile2;
    myfile2.open("susceptibility.txt"); // Columns: L, T, chi

    //initialise temperature range
    iniT = 2.0; finT = 3.0; numT = 101;
    T = linspace(iniT, finT, numT);

    //initialise variables
    double arrayM[dataPoints];
    double chi;
    //double chi_max; //susceptibility at pseudo critical temperature
    //double T_pseudo; //pseudo critical temperature at given L
    //double T_error = (finT-iniT)/(numT-1);
    int L_start = 16; //smallest lattice size
    int L_points = 10; //number of L

    //begin for loop for lattice sizes
    for (int i = 1; i <= L_points; i++) {
        //set lattice size L and recompute N
        L = i*L_start;
        N = L*L;

        F=(int**)malloc(L*sizeof(int*));
        for (int g=0;g<L;g++) F[g]=(int*)malloc(L*sizeof(int));

        //zero chi_max and T_pseudo
        //chi_max = 0.0;
        //T_pseudo = 0.0;

        //begin for loop for temperature range
        for (int j = 0; j < numT; j++) {
            // Loading
           // float pct = 100*i/L_points;
           // float pct2 = 100*j/numT;
            //cout << "\r" << "Susceptibility Loading : " << setprecision(2) << endl << pct << "%, subsection : " << setprecision(2) << pct2 << "%";

            initiStateHot(); // initialise the state
            flipFunction(T[j], nCycles); // perform nCycles of the Monte Carlo sweep (Metropolis)
            tempM = computeMagnetisation();

            // insert the 1st data points after nCycles
            arrayM[0] = abs(tempM);

            //time
            //clock_t tStart = clock();

            //begin for loop of secondary phase to get multiple values of magnetisation for the current temperature
            for (int k = 1; k < dataPoints; k++){
                flipFunction(T[j],mCycles);
                arrayM[k] = fabs(tempM);

            }

            //cout << "dataPoints computed in : "<< (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
            //tStart = clock();
            //compute susceptibility at current T


            //print out chi value at current T in susceptibility text
            myfile2 << L << "," << T[j] << "," ;

            chi = computeSusceptibility(arrayM, dataPoints, T[j]);


            myfile2 << chi << "," << endl;

//            //cout << "chi computed in : "<< (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
//            //tStart = clock();
//            //compare with previous value of chi_max and replace T_pseudo if necessary
//            if (chi > chi_max) {
//                chi_max = chi;
//                T_pseudo = T[j];
//            }
//
       }

       // myfile << L << "," << T_pseudo << "," << T_error << endl;
   }

    //end of sequence, close file and print out that the computations are finished.
    cout << "\r" << "Susceptibility Loading : 100%" << endl;
    if(myfile.is_open()){
        myfile.close();
        cout << "finished susceptibility" << endl;
    }
}

void susceptibilityData2(int L_arg) {
     /*
    Investigating the susceptibility, chi, of the sample over a range of temperatures.
    This will give out data for various lattice sizes L around the critical temperature.
    Since we only really care about the region around T_crit we will set iniT = 2 and finT = 3
    and compute enough points so that fitting a curve will not be necessary.
    We will also print out only the maximum value of chi for the given L.
    Important:
    The function will append the computed data on the files defined below
    */
    //set lattice size L and recompute N
    L = L_arg;
    N = L*L;
    cout << "L = " << L << endl;

    //initialise the files we will be writing in
    //ofstream myfile;
    //myfile.open("Tpseudo_final.txt", ios_base::app); // Columns: L, T_pseudo, T_error

    ofstream myfile2;
    string filename = "susceptibility_"+to_string(L)+".txt";
    filename_rename_if_exists(filename);
    myfile2.open(filename); // Columns: L, T, chi

    //initialise temperature range
    iniT = 2.0; finT = 2.6; numT = 61;
    T = linspace(iniT, finT, numT);

    //initialise variables
    double arrayM[dataPoints];
    double chi; //susceptibility for plot at given L
    //double chi_max; //susceptibility at pseudo critical temperature
    //double T_pseudo; //pseudo critical temperature at given L
    //double T_error = (finT-iniT)/(numT-1);

    double chiError; //chi error bars

    //allocate memory
    delete[] F;
    F=(int**)malloc(L*sizeof(int*));
    for (int g=0;g<L;g++) F[g]=(int*)malloc(L*sizeof(int));

    //zero chi_max and T_pseudo
    //chi_max = 0.0;
    //T_pseudo = 0.0;

    //begin for loop for temperature range
    for (int j = 0; j < numT; j++) {
        //time
        //clock_t tStart = clock();

        //smart calculations
        //if (T[j]<2.20 && (j%10)!=0) {continue;}
        //if (T[j]>2.40 && (j%10)!=0) {continue;}

        // Loading
        float pct = 100*j/numT;
        //cout << "\r" << "Susceptibility Loading : " << setprecision(2) << pct << "%";

        initiStateHot(); // initialise the state
        flipFunction(T[j], nCycles); // perform nCycles of the Monte Carlo sweep (Metropolis)
        tempM = computeMagnetisation();

        // insert the 1st data points after nCycles
        arrayM[0] = fabs(tempM);

        //time
        //cout << "nCycles computed in : "<< (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
        //tStart = clock();

        //begin for loop of secondary phase to get multiple values of magnetisation for the current temperature
        for (int k = 1; k < dataPoints; k++){
            cout << "\r" << "Susceptibility Loading : " << setprecision(2) << pct << "%" << " ( " << k << " / " << dataPoints << " )     ";
            flipFunction(T[j],mCycles);
            arrayM[k] = fabs(tempM);
        }

        //time
        //cout << "dataPoints computed in : "<< (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
        //tStart = clock();

        //compute susceptibility at current T
        chi = computeSusceptibility(arrayM, dataPoints, T[j]);

        //chiError =  chiArray(arrayM, 100, T[j]);
        chiError = bootstrap_error(arrayM, dataPoints, 128)*N/T[j];
        double *autocorr;
        autocorr = autocorrelation(arrayM, dataPoints);
        int decor_t = 0;
        for (int i = 0; i < dataPoints; i++){
           if (autocorr[i] <= exp(-2)){
               decor_t = i;
               break;
           }
        }
        delete[] autocorr;
        chiError *= sqrt(1+decor_t);

        //time
        //cout << "chi computed in : "<< (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
        //tStart = clock();

        //print out chi value at current T in susceptibility text
        myfile2 << T[j] << " " << chi << " " << chiError << endl;
    }

    //myfile << L << "," << T_pseudo << "," << T_error << endl;

    //end of sequence, close file and print out that the computations are finished.
    cout << "\r" << "Susceptibility Loading : 100%" << endl;
    if(myfile2.is_open()){
        myfile2.close();
        cout << "finished susceptibility" << endl;
    }
}

quantity susceptibilityData3(double T_mt) {
    //initialize
    double* arrayM;
    arrayM = new double[dataPoints];
    double chi;
    double chiError;
    double tempM_mt;
    int *F_mt;
    F_mt = new int[V];
    //compute magnetisation
    initiStateHotND_mt(F_mt);
    flipFunctionND_mt(F_mt, T_mt, tempM_mt, nCycles);
    tempM_mt = computeMagnetisationND_mt(F_mt);
    arrayM[0] = fabs(tempM_mt);
    for (int k = 1; k < dataPoints; k++){
        flipFunctionND_mt(F_mt, T_mt, tempM_mt, mCycles);
        arrayM[k] = fabs(tempM_mt);
    }
    //compute susceptibility value and error
    chi = computeSusceptibility(arrayM, dataPoints, T_mt);
    chiError = bootstrap_error(arrayM, dataPoints, 128)*V/T_mt;
    double *autocorr;
    autocorr = autocorrelation(arrayM, dataPoints);
    int decor_t = 0;
    for (int i = 0; i < dataPoints; i++){
       if (autocorr[i] <= exp(-2)){
           decor_t = i;
           break;
       }
    }
    chiError *= sqrt(1+decor_t);
    //add both to return value
    quantity return_value;
    return_value.value = chi;
    return_value.error = chiError;
    // initialize a file to record magnetisation data if required for an investigation
    // ofstream myfile;
    // string filename = to_string(T_mt)+".txt";
    // filename_rename_if_exists(filename);
    // myfile.open(filename);
    // for(int i = 0; i < dataPoints; i++) {
    //    myfile << arrayM[i] << endl;
    // }
    //wrap up
    delete[] F_mt;
    delete[] autocorr;
    cout << T_mt << ": finished" << endl;
    //return
    return return_value;
}

/*********** Testing N dimensional variation of Metropolis algorithm *************/

void testMetropolisND() {
    /*
    //2D
    L = 16;
    dim = 2;
    V = pow(L,dim);
    F_2 = new int[V];
    initiStateHotND();
    print_2D_array(F_2);
    tempM = computeMagnetisationND();
    cout << tempM << endl;
    flipFunctionND(1.1, 1000);
    print_2D_array(F_2);
    cout << tempM << endl;
    */

    //3D
    L=8;
    dim = 3;
    N = L*L;
    V = pow(L,dim);
    F_2 = new int[V];
    initiStateHotND();

    /*
    //print to console
    printArray3D(F_2);
    tempM = computeMagnetisationND();
    cout << tempM << endl;
    flipFunctionND(1.1, 1000);
    printArray3D(F_2);
    cout << tempM << endl;
    */
    //print to file
    ofstream myfile;
    myfile.open("testMetropolis3D.txt");
    printArraytoFile3D(myfile);
}

void testMetropolisND2() {

    L = 16;
    dim = 3;
    T = linspace(iniT,finT,numT);
    N = L*L;
    V = pow(L,dim);
    F_2 = new int[V];

//    for (int i=55; i<numT; i++) {
//        testMetropolisND2_1(T[i],20000,200);
//    }
    testMetropolisND2_1(T[0],10000,100);
}

void testMetropolisND2_1(double T, int cycles, int configs) {

    string filename = to_string(T)+".txt";
    bool exists = file_exists(filename);

    ofstream myfile;
    myfile.open(filename, ios_base::app);

    string filename2; // for configuration file
    string filename3; // for exact configuration of outlier file
    string filename26 = "initial_config_26.txt";

    if (!exists){
        myfile << T << endl;
    }

    for (int j=0; j<configs; j++) {
        initiStateHotND();
        tempM = computeMagnetisationND();
        flipFunctionND(T,cycles);
        // print magnetisation to file1
        myfile << tempM << endl;
        // print configuration to file2 for autoencoder
        filename2 = to_string(T)+"_"+to_string(j)+".txt";
        printArraytoFile3D_2(filename2, T);
        // print more visual configuration of outliers
        if (fabs(tempM) < 0.9) {
            filename3 = to_string(T)+"_"+to_string(j)+"_outlier"+".txt";
            printArray3D_2(F_2, filename3);
        }
    }

    if(myfile.is_open()){
        myfile.close();
        cout << to_string(T)+".txt finished" << endl;
    }

}

/*********** Using N dimensional variation of Metropolis algorithm to produce configurations for the auto-encoder *************/

void produceConfigurations3D(int L_arg, int cycles, int configs) {
    dim =3;
    L = L_arg;
    T = linspace(iniT,finT,numT);
    V = L*L*L;
    F_2 = new int[V];

    for (int i = 0; i<numT; i++){
        produceConfigurations3D_2(T[i], cycles, configs);
    }
}

void produceConfigurations3D_2(double T, int cycles, int configs) {
    // declare filename strings
    string filename1 = to_string(T)+".txt"; // magnetisation data for T temperature
    string filename2; // configuration file name

    ofstream myfile;
    myfile.open(filename1, ios_base::app);
    if (!file_exists(filename1)){
        myfile << T << endl;
    }

    for (int j = 1; j<=configs; j++) {
        initiStateHotND();
        tempM = computeMagnetisationND();
        flipFunctionND(T,cycles);
        // print magnetisation to file1
        myfile << tempM << endl;
        // print configuration to file2 for autoencoder
        filename2 = to_string(T)+"_"+to_string(j)+".txt";
        printArraytoFile3D_2(filename2, T);
    }

    if(myfile.is_open()){
        myfile.close();
        cout << to_string(T)+".txt finished" << endl;
    }
}

/*********** failed attempt at multithreading *************/

quantity magnetisationData_mt(double T_mt, int confs_to_print) {
    //initialise
    double* arrayM;
    arrayM = new double[dataPoints];
    double aveM;
    double errM;
    double tempM_mt;
    int *F_mt;
    F_mt = new int[V];
    //for configurations
    string conf_filename = "conf_"+to_string(T_mt);
    conf_div = dataPoints/confs_to_print;
    conf_counter = 1;
    //compute magnetisation
    initiStateHotND_mt(F_mt);
    flipFunctionND_mt(F_mt, T_mt, tempM_mt, nCycles);
    tempM_mt = computeMagnetisationND_mt(F_mt);
    arrayM[0] = fabs(tempM_mt);
    printConfiguration(conf_filename+"_"+to_string(conf_counter)+".txt");
    conf_counter++;
    for (int k = 1; k < dataPoints; k++){
        flipFunctionND_mt(F_mt, T_mt, tempM_mt, mCycles);
        arrayM[k] = fabs(tempM_mt);
        if (dataPoints%conf_div==0) {
            printConfiguration(conf_filename+"_"+to_string(conf_counter)+".txt");
        }
    }
    //compute averages and error
    aveM = averageArray(arrayM, dataPoints);
    errM = jackknife_new(arrayM, n_b);
    //add both to return value
    quantity return_value;
    return_value.value = aveM;
    return_value.error = errM;
    //wrap up
    delete[] arrayM;
    delete[] F;
    cout << T_mt << ": finished" << endl;
    //return
    return return_value;
}

/*********** print configuration *************/

void printConfiguration(const string &filename) {
    ofstream file;
    filename_rename_if_exists(filename);
    file.open(filename)
    int spin;
    for (int i =0; i<V; i++) {
        if (F_2[i] == 1) {
            spin = 1;
        }
        else {
            spin = 0;
        }
        file << spin << endl;
    }
    if(file.is_open()){
        file.close();
        cout << filename << endl;
    }
}

/*********** main *************/

int main(int argc, char** argv) {
    // initialise the random number sequence in the main.
    srand(seed);
    cout << "Program Seed: " << seed << endl;
    // allocate space for spins array
    // F=(int**)malloc(L*sizeof(int*));
    // for (int i=0;i<L;i++) F[i]=(int*)malloc(L*sizeof(int));
    // delete[] F;

    // producing data for the report
    bool magn = true;
    bool susc = false;
    if (magn) {
        clock_t tStart = clock();
        cout << "running for magnetisation data" << endl;
        L = 16;
        dim = 2;
        Next2Nearest = 0.5;
        V =(int)(pow(L,dim)+0.5);
        nCycles = 10000;        //number of cycles of Monte Carlo primary
        mCycles = 100;           //number of cycles of Monte Carlo secondary
        dataPoints = 1000;      //total data points
        configurations_to_print = 100; // number of configurations to print
        iniT = 1.0; finT = 4.0; numT = 31;
        // iniT = 2.0; finT = 2.1125; numT = 10;
        // iniT = 2.125; finT = 2.2375; numT = 10;
        // iniT = 2.25; finT = 2.3625; numT = 10;
        // iniT = 2.375; finT = 2.4875; numT = 10;
        // iniT = 2.5; finT = 3.4; numT = 10;
        // iniT = 3.5; finT = 4.4; numT = 10;
        // iniT = 4.5; finT = 4.5; numT = 1;
        T = linspace(iniT, finT, numT);
        cout << "L = " << L << endl;
        cout << "Dimensions = " << dim << endl;
        cout << "Thermalisation cycles = " << nCycles << endl;
        cout << "Number of data points = " << dataPoints << endl;
        cout << "Cycles in between data points = " << mCycles << endl;
        cout << "Temperatures to compute are: " << endl;
        print_T_to_console();
        cout << "Starting computations" << endl;
        vector<double> magn;
        vector<double> err_magn;
        quantity temp;
        // compute data
        for (int i = 0; i<numT; i++) {
            temp = magnetisationData_mt(T[i],configurations_to_print);
            magn.push_back(temp.value);
            err_magn.push_back(temp.error);
        }
        // open file
        ofstream myfile;
        string filename = "magnetisation_"+to_string(L)+".txt";
        cout << "file name is " << filename << endl;
        filename_rename_if_exists(filename);
        cout << "creating magnetisation data file with name " + filename << endl;
        myfile.open(filename);
        // write data on open file
        for(int i = 0; i < numT; i++) {
           myfile << T[i] << " " << magn[i] << " " << err_magn[i] << endl;
        }
        // wrap up
        if (myfile.is_open()){
            myfile.close();
        }
        cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
        delete[] T;
    }
    if (susc) {
        clock_t tStart = clock();
        cout << "running for susceptibility data" << endl;
        L = 16;
        dim = 2;
        V = pow(L,dim);
        nCycles = 30000;        //number of cycles of Monte Carlo primary
        mCycles = 100;           //number of cycles of Monte Carlo secondary
        dataPoints = 20000;      //total data points
        // full time
        iniT = 2.0; finT = 2.59; numT = 60;
        // 9 "threads"
        //iniT = 2.0; finT = 2.09; numT = 10; //thread1
        //iniT = 2.1; finT = 2.19; numT = 10; //thread2
        //iniT = 2.2; finT = 2.236; numT = 10; //thread3
        //iniT = 2.24; finT = 2.276; numT = 10; //thread4
        //iniT = 2.28; finT = 2.316; numT = 10; //thread5
        //iniT = 2.32; finT = 2.356; numT = 10; //thread6
        //iniT = 2.36; finT = 2.396; numT = 10; //thread7
        //iniT = 2.4; finT = 2.49; numT = 10; //thread8
        //iniT = 2.5; finT = 2.59; numT = 10; //thread9
        T = linspace(iniT, finT, numT);
        cout << "L = " << L << endl;
        cout << "Dimensions = " << dim << endl;
        cout << "Thermalisation cycles = " << nCycles << endl;
        cout << "Number of data points = " << dataPoints << endl;
        cout << "Cycles in between data points = " << mCycles << endl;
        cout << "Temperatures to compute are: " << endl;
        print_T_to_console();
        cout << "Starting computations" << endl;
        vector<double> susc;
        vector<double> err_susc;
        quantity temp;
        // compute data
        for (int i = 0; i<numT; i++) {
            temp = susceptibilityData3(T[i]);
            susc.push_back(temp.value);
            err_susc.push_back(temp.error);
        }
        // open file
        ofstream myfile;
        string filename = "susceptibility_"+to_string(L)+".txt";
        cout << "file name is " << filename << endl;
        filename_rename_if_exists(filename);
        cout << "creating susceptibility data file with name " + filename << endl;
        myfile.open(filename);
        // write data on open file
        for(int i = 0; i < numT; i++) {
           myfile << T[i] << " " << susc[i] << " " << err_susc[i] << endl;
        }
        // wrap up
        if (myfile.is_open()){
            myfile.close();
        }
        cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
        delete[] T;
    }

    // failed attempts at multithreading
    // run in terminal:
    // g++ -std=c++14 -pthread autocorrelation.cpp bootstrap.cpp general.cpp ising.cpp jacknife.cpp main.cpp Metropolis.cpp wolff.cpp -o
    bool susc_mt = false;
    bool magn_mt = false;
    bool wolff_mt = false;
    if (susc_mt) {
        clock_t tStart = clock();
        cout << "running for susceptibility on multiple threads" << endl;
        L = 20;
        dim = 3;
        V = L*L*L;
        nCycles = 10000;        //number of cycles of Monte Carlo primary
        mCycles = 50;           //number of cycles of Monte Carlo secondary
        dataPoints = 20000;      //total data points
        cout << "L = " << L << endl;
        cout << "Dimensions = " << dim << endl;
        cout << "Thermalisation cycles = " << nCycles << endl;
        cout << "Number of data points = " << dataPoints << endl;
        cout << "Cycles in between data points = " << mCycles << endl;
        vector<future<quantity>> futures;
        vector<double> susc;
        vector<double> err_susc;
        iniT = 4.0; finT = 5.0; numT = 41;
        T = linspace(iniT, finT, numT);

        for (int i =0; i<numT; i++) {
            cout << T[i] << ": creating thread " << endl;
            futures.push_back(async(susceptibilityData3, T[i]));
        }
        cout << "all threads created" << endl;
        for (auto &e : futures){
            quantity temp = e.get();
            susc.push_back(temp.value);
            err_susc.push_back(temp.error);
        }
        // open file
        ofstream myfile;
        string filename = "susceptibility3D_"+to_string(L)+".txt";
        filename_rename_if_exists(filename);
        cout << "creating susceptibility data file with name " + filename << endl;
        myfile.open(filename);
        // write data on open file
        for(int i = 0; i < numT; i++) {
           myfile << T[i] << " " << susc[i] << " " << err_susc[i] << endl;
        }
        // wrap up
        if (myfile.is_open()){
            myfile.close();
        }
        cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
        delete[] T;
    }
    if (magn_mt) {
        clock_t tStart = clock();
        cout << "running for magnetisation data on multiple threads" << endl;
        L = 48;
        dim = 2;
        V = pow(L,dim);
        nCycles = 40000;        //number of cycles of Monte Carlo primary
        mCycles = 100;           //number of cycles of Monte Carlo secondary
        dataPoints = 1000;      //total data points
        cout << "L = " << L << endl;
        cout << "Dimensions = " << dim << endl;
        cout << "Thermalisation cycles = " << nCycles << endl;
        cout << "Number of data points = " << dataPoints << endl;
        cout << "Cycles in between data points = " << mCycles << endl;
        vector<future<quantity>> futures;
        vector<double> magn;
        vector<double> err_magn;
        iniT = 1.0; finT = 5.0; numT = 41;
        T = linspace(iniT, finT, numT);

        for (int i = 0; i<numT; i++) {
            cout << T[i] << ": creating thread " << endl;
            futures.push_back(async(magnetisationData_mt, T[i]));
        }
        cout << "all threads created" << endl;
        for (auto &e : futures){
            quantity temp = e.get();
            magn.push_back(temp.value);
            err_magn.push_back(temp.error);
        }
        // open file
        ofstream myfile;
        string filename = "magnetisation_"+to_string(L)+".txt";
        filename_rename_if_exists(filename);
        cout << "creating magnetisation data file with name " + filename << endl;
        myfile.open(filename);
        // write data on open file
        for(int i = 0; i < numT; i++) {
           myfile << T[i] << " " << magn[i] << " " << err_magn[i] << endl;
        }
        // wrap up
        if (myfile.is_open()){
            myfile.close();
        }
        cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
        delete[] T;
    }
    if (wolff_mt) {
        clock_t tStart = clock();
        cout << "running for wolff susceptibility on multiple threads" << endl;
        vector<future<quantity>> futures;

        int data_points = 1000;
        vector<double> susceptibilities;
        vector<double> err_susceptibilities;

        double *T;

        T = linspace(2.26, 2.4, data_points);
        for(int i = 0; i < data_points; i++) {
           futures.push_back(async(wolff, T[i]));
        }

        for (auto &e : futures){
           quantity temp = e.get();
           susceptibilities.push_back(temp.value);
           err_susceptibilities.push_back(temp.error);
        }

        string filename ("susceptibility" + to_string(L) + ".txt");
        ofstream myfile;
        myfile.open(filename);

        for(int i = 0; i < data_points; i++) {
           myfile << T[i] << "," << susceptibilities[i] << "," << err_susceptibilities[i] << endl;
        }

        if (myfile.is_open()){
           myfile.close();
        }

        cout << "Program ended in " << (double)(clock()-tStart)/CLOCKS_PER_SEC << " seconds" << endl;
    }

    return 0;
}
