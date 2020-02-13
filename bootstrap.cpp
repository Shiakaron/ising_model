#include "header.h"

double* bootstrap_error(double data[], int array_size, int bin_number, bool correlated){
    /*

    The Bootstap Error Analysis:
    1. create bin_number bins of size array_size
    2. randmoly assign - without replacement - data points in the bins from the data[]
    3. For each bin compute an average
    4. Finally compute the average of the averages and its error
    Repeat for bins of deviation where you randomly assign data[]*data[] to get average and error on deviation

    (not in official bootstrap algorithm)
    For correlated data you can find the autocorrelation time and adjust error accordingly.
    
    */
    //bin_number is the number of total arrays of size array_size to be measured for the bootstrap error
    double *O;                          //the average
    O = new double[bin_number]();
    double *O2;                         //the average of squares
    O2 = new double[bin_number]();
    double *dev;                        //the deviation (chi/C)
    dev = new double[bin_number]();

    int x;
    //calculate average for each bin
    for (int i = 0; i < bin_number; i++){
        //assign without replacement
        for (int j = 0; j < array_size; j++){
            x = rand()%array_size;
            O[i] += data[x];
            O2[i] += data[x]*data[x];
        }
        //calculate average for each bin
        O[i] /= array_size;
        O2[i] /= array_size;
        dev[i] = O2[i] - O[i]*O[i];
    }

    //calculate the averages of the bins
    double av_O = 0;
    double av_dev = 0;
    for (int i = 0; i < bin_number; i++){
        av_O += O[i];
        av_dev += dev[i];
    }
    av_O /= bin_number;
    av_dev /= bin_number;

    //calculates the deviations of each bin from the average of all the bins
    double err_O = 0;
    double err_dev = 0;
    for (int i = 0; i < bin_number; i++){
        err_O += (O[i] - av_O)*(O[i] - av_O);
        err_dev += (dev[i] - av_dev)*(dev[i] - av_dev);
    }
    err_O /= bin_number; err_dev /= bin_number;
    err_O = sqrt(err_O); err_dev = sqrt(err_dev);

    //extra feature to get even better estimate on error when the data are correlated
    if (correlated==true) {
        double *autocorr;
        autocorr = autocorrelation(data, array_size);
        int tau_2e = 0;
        for (int i = 0; i < array_size; i++){
           if (autocorr[i] <= exp(-2)){
               tau_2e = i;
               break;
           }
        }
        double term = sqrt(1+tau_2e);
        err_O *= term; err_dev *= term;
    }

    //returns four values, the average of all arrays and its error and the average deviation and its error.
    double *return_value;
    return_value = new double[4];
    return_value[0] = av_O;
    return_value[1] = err_O;
    return_value[2] = av_dev;
    return_value[3] = err_dev;

    delete[] O;
    delete[] O2;
    delete[] dev;

    return return_value;
}
