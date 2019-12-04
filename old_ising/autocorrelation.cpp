#include "header.h"

double rho_calculation(double arr[], int t, int t_max){
    //rho is the single value at point t for the data in array arr[]
    //essentially, the deviation for the values 0 to t_max-t is multiplied by the deviation of the values from t to t_max
    double xav0 = 0, xavt = 0, r = 0;
    int sum_index = t_max - t;

    //this for loop calculates the two averages of the ranges 0 to t_max-t and t to t_max
    for (int t0 = 0; t0 < sum_index; t0++){
        xav0 += arr[t0];
        xavt += arr[t0 + t];
    }
    xav0 /= sum_index;
    xavt /= sum_index;


    //this loop multiplies the individual deviations from the means and returns the rho value at time t
    for (int t0 = 0; t0 < sum_index; t0++){
        r += (arr[t0] - xav0)*(arr[t0+t]-xavt);
    }
    r /= sum_index;
    return r;
}

double* autocorrelation(double arr[], int t_max){
    // rho(t) = <(O(t') - <O>)(O(t'+t)-<O>)>_t' / <(O-<O>)^2>
    // see Autocorrelation chapter in Computational Physics of Konstantinos N. Anagnostopoulos for more info of how this works
    //define rho, the autocorrelation function, size t_max
    double* rho;
    rho = new double[t_max];

    //for all times below t_max, calculate the value rho
    for (int t = 0; t < t_max; t++){
        rho[t] = rho_calculation(arr, t, t_max);
    }

    //normalise, with the first point being equal to 1
    double norm  = 1/rho[0];

    //multiply all other values of the function by the normalisation
    for (int i = 0; i < t_max; i++){
        rho[i] *= norm;
    }

    //returns the autocorrelation function array
    return rho;
}
