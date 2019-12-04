#include "header.h"
/*
ND stands for N dimensions
*/
void initStateColdND() {
    for (int i = 0; i<V; i++) {
        F_2[i] = 1;
    }
}

void initiStateHotND() {
    for (int i = 0; i<V; i++) {
        F_2[i] = 2*(rand()%2)-1;    //return -1 or 1 randomly (50/50)
    }
    //cout << "I am here3\n";
}

double computeMagnetisationND() {
    double m = 0.0;
    for(int i=0;i<V;i++){
            m += F_2[i];
    }
    m /= V;
    return m;

}

double denergyND(int index) {
    // temp will hold the location of the spin to be added
    int temp;
    double de_total = 0;
    double de_nearest = 0;
    double de_next2nearest = 0;
    double de_externalfield = 0;
    //nearest neighbours
    for (int j = 0; j < dim; j++){
        for (int k = -1; k <= 1; k += 2){
            temp = index + k * pow(L,j);
            if (temp < 0) temp += V;
            if (temp >= V) temp -= V;
            de_nearest += F_2[temp];
        }
    }
    de_nearest *= 2*F_2[index];

    //next to nearest neighbours
    if (Next2Nearest != 0) {
        //dim = 1 is a special case
        if (dim == 1) {
            temp = index-2;
            if (temp < 0) temp += V;
            de_next2nearest += F_2[temp];
            temp = index+2;
            if (temp >= V) temp -= V;
            de_next2nearest += F_2[temp];
        }
        else {
            // dim > 1
            for (int j = 0; j < dim-1; j++){
                for (int g = j+1; g < dim; g++) {
                    temp = index + pow(L,j) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index + pow(L,j) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                    temp = index - pow(L,j) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index - pow(L,j) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                }
            }
        }
        de_next2nearest *= 2*Next2Nearest*F_2[index];
    }
    //external field influence
    if (B != 0) {
       de_externalfield += B*F_2[index];
    }
    //total change in energy
    de_total = de_nearest + de_next2nearest + de_externalfield;

    return de_total;
}

double denergyND_periodic(int index){
     // temp will hold the location of the spin to be added
    int temp;
    double de_total = 0;
    double de_nearest = 0;
    double de_next2nearest = 0;
    double de_externalfield = 0;
    //nearest neighbours
    //LEFT
    if (index%L==0) {temp = index + L - 1;}
    else {temp = index - 1;}
    de_nearest += F_2[temp];
    //RIGHT
    if ((index+1)%L==0) {temp = index - L + 1;}
    else {temp = index + 1;}
    de_nearest += F_2[temp];
    //REST
    for (int j=1; j<dim; j++) {
        temp = index - pow(L,j);
        if (temp < 0) temp += V;
        de_nearest += F_2[temp];
        temp = index + pow(L,j);
        if (temp >= V) temp -= V;
        de_nearest += F_2[temp];
    }
    de_nearest *= 2*F_2[index];

    //next to nearest neighbours
    if (Next2Nearest != 0) {
        //cout << "I am here5\n";
        //dim = 1 is a special case
        if (dim == 1) {
            temp = index-2;
            if (temp < 0) temp += V;
            de_next2nearest += F_2[temp];
            temp = index+2;
            if (temp >= V) temp -= V;
            de_next2nearest += F_2[temp];
        }
        else {
            //dim > 1
            //LEFT
            if (index%L==0) {
                for (int g = 1; g < dim; g++) {
                    temp = index + (L-1) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index + (L-1) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                }
            }
            else {
                for (int g = 1; g < dim; g++) {
                    temp = index - 1 + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index - 1 - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                }
            }
            //RIGHT
            if ((index+1)%L==0) {
                for (int g = 1; g < dim; g++) {
                    temp = index - (L-1) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index - (L-1) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                }
            }
            else {
                for (int g = 1; g < dim; g++) {
                    temp = index + 1 + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index + 1 - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                }
            }
            //REST
            for (int j = 1; j < dim-1; j++){
                for (int g = j+1; g < dim; g++) {
                    temp = index + pow(L,j) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index + pow(L,j) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                    temp = index - pow(L,j) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += F_2[temp];
                    temp = index - pow(L,j) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += F_2[temp];
                }
            }
        }
        de_next2nearest *= 2*Next2Nearest*F_2[index];
    }
    if (B != 0) {
       de_externalfield += B*F_2[index];
    }
    //total change in energy
    de_total = de_nearest + de_next2nearest + de_externalfield;

    return de_total;
}

void flipFunctionND(double T, int cycles) {
    double de = 0;
    //boltzman factor
    double beta = 1.0/T;
    //initialise random number variable
    double x;
    double flip = 0.0;
    int index;
    //for testing outliers
    //int flip_values[cycles];
    //cout << "I am here4\n";
    for(int r =0; r < cycles; r++){

        for(int k = 0; k<V; k++){ //begin attempting to flip

            // randomly select a site
            index = rand()%V;

            //measure energy difference before flip
            //de = denergyND(index);
            de = denergyND_periodic(index);

            if (de <= 0){
                F_2[index] *= -1;
                flip += F_2[index];

            }
            else {
                x = ((double)rand()/RAND_MAX);
                if (x < exp(-beta*de)){
                    F_2[index] *= -1;
                    flip += F_2[index];
                }
            }
        }

        //for testing outliers
        //flip_values[r] = (int)flip;
    }
    tempM += 2*flip/V;

    /*
    //for testing outliers
    if (fabs(tempM)<0.9) {
        //name the file
        int t = 2;
        string filename = to_string(tempM)+"_flip_values.txt";
        bool exists = file_exists(filename);
        while (exists) {
            filename = to_string(tempM)+"_flip_values("+to_string(t)+").txt";
            exists = file_exists(filename);
            t += 1;
        }
        //open the file
        ofstream myfile;
        myfile.open(filename);
        //write the values
        for (int i =0; i<cycles; i++){
            myfile << flip_values[i] << endl;
        }
        //close the file
        if(myfile.is_open()){
            myfile.close();
            cout << filename << endl;
        }
    }
    */
}

/** recreatuing the functions that take array as an argument in order to implement multithreading **/

void initiStateHotND_mt(int arr[]) {
    for (int i = 0; i<V; i++) {
        arr[i] = 2*(rand()%2)-1;    //return -1 or 1 randomly (50/50)
    }
}

double computeMagnetisationND_mt(int arr[]) {
    double m = 0.0;
    for(int i=0;i<V;i++){
        m += arr[i];
    }
    m /= V;
    return m;
}

double denergyND_periodic_mt(int arr[], int index){
     // temp will hold the location of the spin to be added
    int temp;
    double de_total = 0;
    double de_nearest = 0;
    double de_next2nearest = 0;
    double de_externalfield = 0;
    //nearest neighbours
    //LEFT
    if (index%L==0) {temp = index + L - 1;}
    else {temp = index - 1;}
    de_nearest += arr[temp];
    //RIGHT
    if ((index+1)%L==0) {temp = index - L + 1;}
    else {temp = index + 1;}
    de_nearest += arr[temp];
    //REST
    for (int j=1; j<dim; j++) {
        temp = index - pow(L,j);
        if (temp < 0) temp += V;
        de_nearest += arr[temp];
        temp = index + pow(L,j);
        if (temp >= V) temp -= V;
        de_nearest += arr[temp];
    }
    de_nearest *= 2*arr[index];

    //next to nearest neighbours
    if (Next2Nearest != 0) {
        //dim = 1 is a special case
        if (dim == 1) {
            temp = index-2;
            if (temp < 0) temp += V;
            de_next2nearest += arr[temp];
            temp = index+2;
            if (temp >= V) temp -= V;
            de_next2nearest += arr[temp];
        }
        else {
            //dim > 1
            //LEFT
            if (index%L==0) {
                for (int g = 1; g < dim; g++) {
                    temp = index + (L-1) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += arr[temp];
                    temp = index + (L-1) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += arr[temp];
                }
            }
            else {
                for (int g = 1; g < dim; g++) {
                    temp = index - 1 + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += arr[temp];
                    temp = index - 1 - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += arr[temp];
                }
            }
            //RIGHT
            if ((index+1)%L==0) {
                for (int g = 1; g < dim; g++) {
                    temp = index - (L-1) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += arr[temp];
                    temp = index - (L-1) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += arr[temp];
                }
            }
            else {
                for (int g = 1; g < dim; g++) {
                    temp = index + 1 + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += arr[temp];
                    temp = index + 1 - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += arr[temp];
                }
            }
            //REST
            for (int j = 1; j < dim-1; j++){
                for (int g = j+1; g < dim; g++) {
                    temp = index + pow(L,j) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += arr[temp];
                    temp = index + pow(L,j) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += arr[temp];
                    temp = index - pow(L,j) + pow(L,g);
                    if (temp >= V) temp -= V;
                    de_next2nearest += arr[temp];
                    temp = index - pow(L,j) - pow(L,g);
                    if (temp < 0) temp += V;
                    de_next2nearest += arr[temp];
                }
            }
        }
        de_next2nearest *= 2*Next2Nearest*arr[index];
    }
    if (B != 0) {
       de_externalfield += B*arr[index];
    }
    //total change in energy
    de_total = de_nearest + de_next2nearest + de_externalfield;

    return de_total;
}

void flipFunctionND_mt(int arr[], double T_mt, double& tempM_mt, int cycles) {
    double de = 0;
    //boltzman factor
    double beta = 1.0/T_mt;
    //initialise random number variable
    double x;
    double flip = 0.0;
    int index;


    for(int r =0; r < cycles; r++){
        //begin attempting to flip
        for(int k = 0; k<V; k++){

            // randomly select a site
            index = rand()%V;

            de = denergyND_periodic_mt(arr, index);

            if (de <= 0){
                arr[index] *= -1;
                flip += arr[index];

            }
            else {
                x = ((double)rand()/RAND_MAX);
                if (x < exp(-beta*de)){
                    arr[index] *= -1;
                    flip += arr[index];
                }
            }
        }
    }
    tempM_mt += 2*flip/V;

}
