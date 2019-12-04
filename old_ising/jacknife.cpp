#include "header.h"


double jackknife(int n_b, double arr[]) {
    /*
     function that returns the jackknife error when given the number of bins and the sample array
     */
    double sigmaJack=0.0; // initialise return variable
    int n = dataPoints; // just to make things clearer
    int b = n-n/n_b; // initialise bin size
    double avgO[n_b]; // initialise array holding the average value for each bin
    double avg_avgO; // initialise the variable name for the average value of the averages of bins :p
    double tempO[b]; // initialise recyclable array holding the values of a bin
    
    //create a for loop that updates the i value of avgO on each iteration
    for (int i=0; i<n_b; i++) {
        
        // evaluate jump points
        int j_1 = i*n/n_b;
        int j_2 = (i+1)*n/n_b;
        
        // 1st for loop copies data up to j_1
        for (int j=0; j<j_1; j++) {
            tempO[j] = arr[j];
        }
        // jump point
        // 2nd for loop copies data from j_2
        for (int j=j_2; j<n; j++) {
            tempO[j-n/n_b] = arr[j];
        }
        avgO[i] = averageArray(tempO, b);
    }
    //compute the average of averages :p
    avg_avgO = averageArray(avgO, n_b);
    
    //create a for loop that computes the sum (=sigmaJack^2)
    double sum = 0.0;
    for (int i = 0; i<n_b; i++) {
        sum += (pow((avgO[i]-avg_avgO),2));
    }
    
    sigmaJack = sqrt(fabs(sum));
    return sigmaJack;
}

double* jackknife2(int n_b, double arr[]) {
    /*
     function that returns a pointer to a size 2 array that includes
     at position [0] the average of the sample
     at position [1] the jackknife error
     both computed in the same way as in the jackknife() function
     */
    double *ret;
    ret = new double[2];
    int n = dataPoints; // just to make things clearer
    int b = n-n/n_b; // initialise bin size
    double avgO[n_b]; // initialise array holding the average value for each bin
    double tempO[b]; // initialise recyclable array holding the values of a bin
    
    //create a for loop that updates the i value of avgO on each iteration
    for (int i=0; i<n_b; i++) {
        
        // evaluate jump points
        int j_1 = i*n/n_b;
        int j_2 = (i+1)*n/n_b;
        
        // 1st for loop copies data up to j_1
        for (int j=0; j<j_1; j++) {
            tempO[j] = arr[j];
        }
        // jump point
        // 2nd for loop copies data from j_2
        for (int j=j_2; j<n; j++) {
            tempO[j-n/n_b] = arr[j];
        }
        avgO[i] = averageArray(tempO, b);
    }
    //compute the average of averages :p
    ret[0] = averageArray(avgO, n_b);
    
    //create a for loop that computes the sum (=sigmaJack^2)
    double sum = 0.0;
    for (int i = 0; i<n_b; i++) {
        sum += (pow((avgO[i]-ret[0]),2));
    }
    
    ret[1] = sqrt(fabs(sum));
    
    return ret;
}

void jackknifeErrorAnalysis() {
    /*
     jackknife error analysis to decide number of bins to be used (8)
     Computes data to plot error vs number of bins
     */
    // initialise temperature array
    T = linspace(iniT, finT, numT);
    //jackknife error analysis
    ofstream myfile2;
    myfile2.open("jackknifeM.txt");
    
    ofstream myfile3;
    myfile3.open("jackknifeE.txt");
    
    double arrayM[dataPoints]; // initialise reusable arrays to insert data points, these arrays are used to compute an average and a sigma at each Temperature.
    double arrayE[dataPoints]; // the +1 is for the 1st data point after nCycles and the rest is for the secondary phase with mCycles in between them
    
    double sigmaJackM;
    double sigmaJackE;
    
    // begin for-loop to compute values for each Temperature.
    for (int i =0; i<numT; i++) {
        // Loading
        float pct = 100*i/numT;
        cout << "\r" << "Jackknife Loading : " << setprecision(2) << pct << "%";
        
        myfile2 << "T = " << T[i] << "\n";
        myfile3 << "T = " << T[i] << "\n";
        
        initiStateHot(); // initialise the state
        flipFunction(T[i], nCycles); // perform nCycles of the Monte Carlo sweep (Metropolis)
        
        // insert the 1st data points after nCycles
        arrayM[0] = computeMagnetisation();
        arrayE[0] = computeEnergy();
        
        
        // for loop for secondary phase to get an average for Energy and Magnetisation
        for (int j = 1; j< dataPoints; j++){
            flipFunction(T[i],mCycles);
            arrayM[j]=computeMagnetisation();
            arrayE[j]=computeEnergy();
        }
        
        for (int j=dataPoints; j>1; j--) {
            // j = number of bins
            if ((dataPoints)%j!=0) {continue;}
            
            sigmaJackM = jackknife( j, arrayM);
            sigmaJackE = jackknife( j, arrayE);
            
            // bin numbers, sigmajack
            myfile2 << j << ",";
            myfile2 << sigmaJackM << "\n";
            myfile3 << j << ",";
            myfile3 << sigmaJackE << "\n";
        }
        myfile2 << "\n";
        myfile3 << "\n";
    }
    //============================================================
    //end of sequence, close file and print out that the computations are finished.
    cout << "\r" << "Jackknife Loading : 100%" << endl;
    
    if(myfile2.is_open()){
        myfile2.close();
        cout << "finished sigmaJackM" << endl;
    }
    
    if(myfile3.is_open()){
        myfile3.close();
        cout << "finished sigmaJackE" << endl;
    }
    
}

void jackknifeErrorAnalysis2() {
    /*
     Another jackknife error analysis to decide number of bins to be used (8).
     This time we plotted average Magnetisation vs Bins with their error bars.
     */
    //initialise temperature
    T = linspace(1.2, 3.0, 7);
    
    //magnetisation as a function of bins (?)
    ofstream myfile4;
    myfile4.open("MagVsBins.txt");
    
    double arrayM[dataPoints];
    double *values; // average and error
    
    for (int i = 0; i<7; i++) {
        // Loading
        float pct = 100*i/7;
        cout << "\r" << "Loading : " << setprecision(2) << pct << "%";
        
        initiStateHot();
        flipFunction(T[i], nCycles);
        
        arrayM[0] = computeMagnetisation();
        
        for (int j = 1; j<dataPoints; j++){
            flipFunction(T[i],mCycles);
            arrayM[j]=computeMagnetisation();
        }
        
        for (int j=dataPoints; j>1; j--) {
            // j = number of bins
            if ((dataPoints)%j!=0) {continue;}
            
            values = jackknife2(j, arrayM);
            
            // bin numbers, sigmajack, aveJack
            myfile4 << j << "," << values[0] << "," << values[1] << "\n";
        }
    }
    //end of sequence, close file and print out that the computations are finished.
    cout << "\r" << "Jackknife (2) Loading : 100%" << endl;
    if(myfile4.is_open()){
        myfile4.close();
        cout << "finished MagVsBins" << endl;
    }
}

double jackknife_new(double data[], int n_b){
    //retyrbs error if division does not return an integer
    if (dataPoints%n_b != 0)
        cout << "Error: data points not divisible by bin number" << endl;
    
    int n_removed = dataPoints/n_b;                     //number of data removed from each bin
    double sum = sumArray(data, dataPoints);            //sum of all the daata
    double sums_array[n_b];
    
    //set all array values to the sum
    for (int i = 0; i < n_b; i++){
        sums_array[i] = sum;
    }
    
    //subtract the ith data point from the desired set's sum
    for (int i = 0; i < dataPoints; i++){
        sums_array[(int)floor(i/n_removed)] -= data[i];
    }
    
    //average each set
    for (int i = 0; i < n_b; i++){
        sums_array[i] /= (dataPoints-n_removed);
    }
    
    //find the average of all sets
    double array_average = averageArray(sums_array, n_b);
    sum = 0.0; //reduce reuse recycle
    
    //calculate the deviation of each set from the mean of all sets
    for (int i = 0; i<n_b; i++) {
        sum += (pow((sums_array[i]-array_average),2));
    }
    
    //return jacknife error
    return sqrt(fabs(sum));
}
