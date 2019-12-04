#include "header.h"

double* linspace(double a, double b, int N){
    // similar to linspace of Matlab. Help to create an array of equally spaced N elements from a to b.
    double D =(b-a)/(N-1); // spacing
    double *ret;           // pointer
    ret = new double[N];   // initialise array

    //set in between
    for (int i = 0; i<N; i++) {
        ret[i] = a+D*i;
    }
    ret[N-1] = b;

    return ret;
}

double sumArray(double arr[], int siz) {
    // returns the sum of the elements of the array of size siz.
    double sum = 0.0;
    for (int i=0; i<siz; i++) {
        sum += arr[i];
    }
    return sum;
}

double averageArray(double arr[], int siz) {
    // returns the average of the elements of the array of size siz.
    double sum = sumArray(arr, siz);
    double average = sum/siz;
    return average;
}

double calcSigma(double arr[], int siz){
    //returns the standard deviation of the elements of the array of size siz.
    double sigma = 0.0;
    double sum = 0.0;
    double ave;

    // average of m and e
    //array holding each value
    ave = averageArray(arr,siz);

    for(int i =0; i <  siz; i++){
        sum += pow((arr[i] - ave),2);
    }

    sigma = sqrt(sum)/(siz);

    return sigma;
}

void printArray() {
    // prints the configuration of the spins on the console in a matrix form
    for (int j =0; j<L; j++){
        for(int i=0;i<L;i++){
            cout << (F[i][j]+1) << " ";
        }
        cout << endl;
    }
}

void printArraytoFile(ofstream &file) {
    // prints the configuration of the spins on the file in a matrix form
    for (int j =0; j<L; j++){
        for(int i=0;i<L;i++){
            file << (F[i][j]+1) << ",";
        }
        file << endl;
    }
    file << endl;
}

void printArraytoFile2(ofstream &file) {
    // prints the configuration of the spins on the file in column form.
    // x,y,up/down.
    // x and y start from 1
    for (int j =0; j<L; j++){
        for(int i=0;i<L;i++){
                file << (i+1) << " " << (j+1) << " " << F[i][j] << endl;
        }
    }
}

void printArray3D(int arr[]) {
    if (dim != 3 ) {
        cout << "The dimensions need to be 3 for the function printArray3D() to work as intended" << endl;
    }
    else {
        for (int i = 0; i < V; i++){
            if (i % N == 0)
                cout << "layer " << i/N << endl;
            if (arr[i] == 1)
                cout << "1" << " ";
            else if (arr[i] == -1)
                cout << "0" << " ";
            if ((i+1) % L == 0)
                cout << endl;

        }
    }
}

void printArray3D_2(int arr[], const string &filename) {
    if (dim != 3 ) {
        cout << "The dimensions need to be 3 for the function printArray3D() to work as intended" << endl;
    }
    else {
        ofstream file;
        file.open(filename);

        for (int i = 0; i < V; i++){
            if (i % N == 0)
                file << "layer " << i/N << endl;
            if (arr[i] == 1)
                file << "1" << " ";
            else if (arr[i] == -1)
                file << "0" << " ";
            if ((i+1) % L == 0)
                file << endl;
        }
        if(file.is_open()){
            file.close();
            cout << filename << endl;
        }
    }
}

void printArraytoFile3D(ofstream &file) {
    if (dim != 3 ) {
        cout << "The dimensions need to be 3 for the function printArraytoFile3D() to work as intended" << endl;
    }
    else {
        // x, y, z, up/down
        int spin;
        for (int i =0; i < V; i++) {
            if (F_2[i] == 1) {
                spin = 1;
            }
            else {
                spin = 0;
            }
            int x = i%L;
            int y = floor(i/L);
            y = y%L;
            int z = floor(i/N);
            /*
            add +1 to x y and z if you want them to start from 1
            */
            file << x << " " << y << " " << z << " " << spin << endl;
        }
    }
}

void printArraytoFile3D_2(const string &filename, double T) {
    ofstream file;
    file.open(filename);
    //print temperature on first line
    file << T << endl;
    //proceed with spins
    int spin;
    for (int i =0; i < V; i++) {
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

bool file_exists(const string& name) {
    ifstream infile(name);
    return infile.good();
}

void filename_rename_if_exists(string& filename) {
    bool exists = file_exists(filename);
    // find dot index
    if (exists) {
        string name;
        string ending;
        int index = filename.length()-1;
        bool found_dot = false;
        while (!found_dot){
            if (filename[index] == '.') {
                found_dot = true;
            }else if(index == 0){
                cout << "no dot \".\" found" << endl;
            }else {
                index--;
            }
        }
        int dot_index = index;
        // break the string apart
        name = filename.substr(0,dot_index);
        //cout << name << endl;
        ending = filename.substr(dot_index);
        //cout << ending << endl;
        // add the (number) in between the two parts
        int t = 0;
        while (exists) {
            string number = "("+to_string(t)+")";
            filename = name+number+ending;
            cout << "file name changed to " << filename << endl;
            exists = file_exists(filename);
            t++;
        }
        //cout << filename << endl;
    }

}

void print_T_to_console() {
    cout << "[ ";
    for (int i = 0; i<numT-1; i++) {
        cout << T[i] << ", ";
    }
    cout << T[numT-1] << " ]" << endl;
}
