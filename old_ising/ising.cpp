#include "header.h"

void initiStateCold() {
    // initialise the configuration with all spins aligned
    for (int i =0; i<L; i++){
        for(int j=0;j<L;j++){
            F[i][j] = 1;   //set all spins aligned
        }
    }
}

void initiStateHot(){
    // initialise the configuration with random spins
    int randNum=0;
    int i,j;                           //identify node in function Q
    for (i =0; i<L; i++){
        for(j=0;j<L;j++){
            randNum = rand()%2;       //generate num; either 0 or 1
            F[i][j] = 2*randNum-1;    //return -1 or 1 randomly (50/50)
        }
    }
}

double computeEnergy() {
    // compute the energy of the configuration
    double e = 0.0;
    for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
             //compute energy
            e -= (F[i][j]) * ((F[i][(j+1)%L] + F[(i+1)%L][j] + Next2Nearest*(F[i?(i-1):(L-1)][(j+1)%L] + F[i?(i-1):(L-1)][j?(j-1):(L-1)])));
        }
    }

    //normalize
    e /= (2*N);

    return e;
}

double denergy(int i, int j){
    //using markov process and metropolis algorithm energy difference is sum
    //using periodic boundaries
    double de = 2*(F[i][j])*(F[i?(i-1):(L-1)][j]+F[(i+1)%L][j]+F[i][j?(j-1):(L-1)]+F[i][(j+1)%L]);
    //helical boundaries
    // int x = 2*(F[i][j])*(F[i?(i-1):(L-1)][i?j:(j+4)%5]+F[(i+1)%L][i==L-1?(j+1)%L:j]+F[i][j?(j-1):(L-1)]+F[i][(j+1)%L]);
    de += 2*Next2Nearest*(F[i?(i-1):(L-1)][(j+1)%L] + F[i?(i-1):(L-1)][j?(j-1):(L-1)]+F[(i+1)%L][(j+1)%L] + F[(i+1)%L][j?(j-1):(L-1)]);
    de += B*F[i][j];


    return de;
    //use periodic for small lattices and helical for large ones for faster computing
}

void flipFunction(double T, int cycles) {
    /*
    A central function to the program that encompasses the Metropolis algorithm.
    A random spin is selected in the configuration and the function considers the energy change in flipping the spin.
    If the energy change is negative then the spin is flipped. If it's positive then there's a probability for it
    to flip according to a boltzman factor.
    The function takes as arguments the temperature the system is evolving in and the number of cycles it will run for.
    For cycles=1, N random points are considered for flipping where N are the number of spins in the square system.
    Another important part of the function is that it updates the magnetisation of the system when it finishes its cycles
    by using a temporary variable "flip" which it increments or decrements according to the spin that it flips. ie if
    the spin goes from +1 to -1 then the flip function is decremented by 1 and the magnetisation global variable tempM
    is changed by 2*flip, or in this case, -2.
    Also, in the case of the Ising Model we have only 5 values for the energy changes, {-8,-4,0,4,8}, so it would be more
    efficient to save two variables for the probabilities of flipping in the case of 4 and 8 instead of calculating the
    boltzman factor exp(-beta*de) every time. We used to have a flipFunction() like that but it must have been deleted somewhere
    along the way leaving behind the more generic case which still works just fine but it must be slower. The Ising specific
    flipFunction used the switch (case) function of C++ language.
    */
    double de = 0;
    //boltzman factor
    double beta = 1.0/T;
    double x;
    double flip = 0.0;
    int i, j;

    for(int r =0; r < cycles; r++){
        //initialise random number variable

        for(int k = 0; k<N; k++){ //begin attempting to flip

            // randomly select a site
            i = rand()%L;
            j = rand()%L;

            de = denergy(i,j); //measure energy difference before flip

            if (de <= 0){
                F[i][j] *= -1;
                flip += F[i][j];

            }
            else if (de > 0){
                x = ((double)rand()/RAND_MAX);
                if (x < exp(-beta*de)){
                    F[i][j] *= -1;
                    flip += F[i][j];
                }
            }
        }
    }
    tempM = (tempM + (2*flip/N));
}

double computeMagnetisation() {
    /*
    this function considers all the spins of the current configuration and adds them up to get the
    total magnetisation, then divides by N to get the average magnetisation. This function is called
    once at the start of the program to initialise the magnetisation global variable tempM which is later
    updated by the flipFunction.
    */
    double m = 0.0;
    for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
            m += F[i][j]; //magnetization; sum of spins
        }
    }

    //normalize
    m /= N;
    //m = fabs(m);

    return m;
}

double computeSusceptibility(double arr[], int siz, double T) {
    /*
    chi = beta*N*(<m^2>-<m>^2)
    */
    //boltzman factor
    double beta = 1/T;
    // < m >
    double aveM = averageArray(arr, siz);
    // < m^2 >
    double aveM_sq;
    double sum = 0.0;

    for (int i = 0; i<siz; i++) {
        sum += arr[i]*arr[i];
    }
    aveM_sq=sum/siz;

    return beta*N*(aveM_sq-aveM*aveM);
}

double chiArray(double arrayM[], int binSize, double T){
    /*
    This function is used to compute the error in the susceptibility calculation. The idea is to split the magnetisation
    calculations in arrays of size binSize, compute averages from these arrays and then calculate the error of these averages.
    The arrays are created by taking adjacent data in the arrayM, thus a better function can be created where you use the mod
    function to get more spaced out data. eg for 1000 datapoints in arrayM and binSize 50 you can create 20 bins with every
    20th element of arrayM. This is will be a bit more expensive to compute since you need to traverse the whole array for each bin instead
    of stopping after 50 elements but the correlation of the datapoints in each bin will decrease dramatically.
    */
    double tempM[binSize]; //temp store values
    //break my array into smaller arrays
    int n = dataPoints/binSize;
    double chi[n]; //store values of chi
    int k;
    for(int i = 0; i < n; i++ ){
        k=0;
        for(int j=i; j < dataPoints; j +=n ){
            tempM[k] = arrayM[j];
            k++;
        }
        chi[i] = computeSusceptibility(tempM, binSize, T);
     }

    double error = calcSigma(chi, n);

    return error;
}
