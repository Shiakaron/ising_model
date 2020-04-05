#include "header.h"

double* linspace(double a, double b, int N){
    // similar to linspace of Matlab. Help to create an array of equally spaced N elements from a to b.
    // returns POINTER to that array
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

bool file_exists(const string& p_file) {
    // 1. Opens the file identified by argument "name", associating it with the stream object,
    // 2. Check whether state of stream is good (public member function )
    ifstream infile(p_file);
    return infile.good();
}

void filename_rename_if_exists(string& filename, string& folder) {
    // check if a file with name=filename already exists
    bool exists = file_exists(folder+"\\"+filename);
    // if it exists, we rename the file to "filename(number).filetype"
    if (exists) {
        // initialise useful variables
        string name; // filename
        string ending; // .filetype
        int index = filename.length()-1; // initialise index to the last element
        bool found_dot = false;
        // find dot index
        while (!found_dot) {
            if (filename[index] == '.') {
                found_dot = true; // the dot's index is found => break from while loop
            }else if(index == 0) {
                cout << "no dot \".\" found" << endl; // not treated as an error but this is never expected to occur for .txt files
            }else {
                index--; // dot's index not found => repeat while loop
            }
        }
        int dot_index = index; // set dot's index
        // break the string apart at the dot's index
        name = filename.substr(0,dot_index);
        ending = filename.substr(dot_index);
        // add the (number) in between the two parts
        int t = 1;
        while (exists) {
            string number = "("+to_string(t)+")";
            filename = name+number+ending;
            exists = file_exists(folder+"\\"+filename); // if the file still exists increment t
            t++;
        }
        cout << "file name changed to: " << filename << endl;
    }
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

observable compute_average_and_sigma(double arr[], int siz) {
    observable O;
    //average first
    O.value = averageArray(arr,siz);
    //sigma
    double sum = 0.0;
    for(int i=0; i<siz; i++){
        sum += (arr[i] - O.value)*(arr[i] - O.value);
    }
    O.error = sqrt(sum)/(siz);
    return O;
}

void print_T(int siz) {
    // print on the console the Temperature array (T) in a readable format. Takes as argument the size of the array
    cout << "Temperatures to compute are:" << endl;
    cout << "[ ";
    for (int i = 0; i<siz-1; i++) {
        cout << T[i] << ", ";
    }
    cout << T[siz-1] << " ]" << endl;
}

void print_mapOfNearest(int max){
    // print on the console the nearest neigbours. For large systems limit the output to max number of sites.
    cout << "Map of Nearest:" << endl;
    vector<int> temp_vec;
    for (int i = 0; i<N; i++) {
        cout << i << ": ";
        temp_vec = mapOfNearest[i];
        for (int j = 0; j<temp_vec.size(); j++) {
            cout << temp_vec[j] << " ";
        }
        cout << endl;
        // skip some elements in the map using the max argument
        if (i == max/2) {
            cout << "..." << endl;
            i = N - max/2;
        }
    }
}

void print_mapOfNext2Nearest(int max){
    // print on the console the next-to-nearest neigbours. For large systems limit the output to max number of sites.
    cout << "Map of Next2Nearest:" << endl;
    vector<int> temp_vec;
    for (int i = 0; i<N; i++) {
        cout << i << ": ";
        temp_vec = mapOfNext2Nearest[i];
        for (int j = 0; j<temp_vec.size(); j++) {
            cout << temp_vec[j] << " ";
        }
        cout << endl;
        // skip some elements in the map using the max argument
        if (i == max/2) {
            cout << "..." << endl;
            i = N - max/2;
        }
    }
}

void print_all_parameters(int thermalisationCycles, int dataPoints, int spacingCycles, int numT, double Temp) {
    // print on the console all the parameters
    cout << "Program Seed: " << seed << endl;
    cout << "Dimensions: " << dim << endl;
    cout << "Lattice size: " << L << endl;
    cout << "Next-to-nearest neighbour interaction: " << n2n << endl;
    cout << "External Field: " << H << endl;
    cout << "Thermalisation cycles: " << thermalisationCycles << endl;
    cout << "Number of data points: " << dataPoints << endl;
    cout << "Cycles in between data points: " << spacingCycles << endl;
    if (numT > 0) {
        print_T(numT);
    }
    else {
        cout << "Temperature of system: " << Temp << endl;
    }
}

void print_spins() {
    // print on the console the spin for each site
    cout << "Spins:" << endl;
    cout << "index" << "\t" << "spin" << endl;
    for (int i=0; i<N; i++) {
        cout << i << "\t" << spins[i] << endl;
    }
}

void print_array(double arr[], int siz) {
    // print on the console an array of doubles and size siz.
    cout << "[ ";
    for (int i = 0; i<siz-1; i++) {
        cout << arr[i] << ", ";
    }
    cout << arr[siz-1] << " ]" << endl;
}

void print_spins_2D() {
    // print on the console the spin array with the following index format (L = 3)
    // 6 7 8
    // 3 4 5
    // 0 1 2
    // this was used to check the code and is only useful for 2D systems
    if (dim == 2) {
        int s;
        cout << "Spins: " << endl;
        for (int j=L-1; j>=0; j--) {
            for (int k=0; k<L; k++) {
                s = spins[j*L+k];
                if (s == 1) {
                    cout << 1 << " ";
                }
                else {
                    cout << 0 << " ";
                }
            }
            cout << endl;
        }
    }
    else {
        cout << "The function print_spins_2D() only works properly with a 2D system.\n";
        cout << "Currently dim = " << dim << endl;
    }
}

int user_integer_input(int min, int max) {
    /*
    This function performs the necessary checks on the user's input. A piece of this code was found online and was added upon check for all invalid inputs. Another addition was made to limit the inputs between integers min and max.
    */
    int x; // input variable
    bool valid = false;
    while (!valid) // loop continually until valid input received
    {
        if (! (cin >> x) ) {  // check stream state
            // eof() or bad() break read loop
            if (cin.eof() || cin.bad()) {
                cerr << "(user canceled or unreconverable error)\n";
                return 1;
            }
            else if (cin.fail()) { // if failbit
                cerr << "error: invalid input.\n";
                cin.clear(); // clear failbit
                // extract any characters that remain unread
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
            }
        }
        else {  // on succesful read of int
            // check if there are more characters
            int c = cin.peek();
            // cout << "cin.peek= " << c << endl;
            if ( c == EOF ) {
                cerr << "(user canceled or unreconverable error)\n";
                return 1;
            }
            else if (c!=10 && c!=32) {
                // 10 == no character after integer input
                // 32 == space after integer input
                cerr << "error: invalid input.\n";
                cin.clear(); // clear failbit
                // extract any characters that remain unread
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
            }
            // check if int is in range min-max
            else if ( x < min || x > max) {
                cerr << "error: input not in range.\n";
                cin.clear(); // clear failbit
                // extract any characters that remain unread
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
            }
            else {
                valid = true;  // then break read loop
            }
        }
    }
    return x;
}
