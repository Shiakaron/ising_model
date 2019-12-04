#include "header.h"

double bootstrap_error(double data[], int array_size, int bin_number){
    //bin number is the number of total arrays of size array_size to be measured for the bootstrap error
    double *O;                          //the average
    O = new double[bin_number]();
    double *O2;                         //the average of squares
    O2 = new double[bin_number]();
    double *dev;                        //the deviation (or chi)
    dev = new double[bin_number]();

    int x = 0;
    //calculate average for each bin
    for (int i = 0; i < bin_number; i++){
        for (int j = 0; j < array_size; j++){
            x = rand()%array_size;
            O[i] += data[x];
            O2[i] += data[x]*data[x];
        }
        O[i] /= array_size;
        O2[i] /= array_size;
        dev[i] = O2[i] - O[i]*O[i];
    }

    //calculate the averages of the bins
    double av_mean = 0;
    double av_dev = 0;
    for (int i = 0; i < bin_number; i++){
        av_mean += O[i];
        av_dev += dev[i];
    }
    av_mean /= bin_number;
    av_dev /= bin_number;

    //calculates the deviations of each bin from the average of all the bins
    double err_mean = 0;
    double err_dev = 0;

    for (int i = 0; i < bin_number; i++){
        err_mean += pow((O[i] - av_mean),2);
        err_dev += pow((dev[i] - av_dev),2);
    }
    err_mean /= bin_number;
    err_dev /= bin_number;

    err_mean = sqrt(err_mean);
    err_dev = sqrt(err_dev);

    //returns four values, the average of all arrays and its error and the average deviation and its error.
    delete[] O;
    delete[] O2;
    delete[] dev;
    return err_dev;

}
