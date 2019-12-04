#include "header.h"

//CYCLIC USES THESE CONDITIONS:
/* NEED TO CHECK IF WE ON A BOUNDARY:
 UPPER BOUNDARY: floor(i/L) = 0
 LOWER BOUNDARY: ceil((i+1)/L) = L
 LEFT BOUNDARY: i % L = 0
 RIGHT BOUNDARY: (i+1) % L = 0
 */

double average_grain_size(int S[]){
    int *flag;
    flag = new int[N];

    int index = 0;
    int cluster_size = 0;

    vector<int> stack;
    vector<int> sizes;

    for (int i = 0; i < N; i++){
        if (flag[i] == 1) continue;
        else{
            flag[i] = 1;
            stack.push_back(i);
            cluster_size = 1;

            while (stack.size()>0){
                index = stack.back();
                stack.pop_back();

                int up = index - L;
                if (up < 0) up += N;

                int down = index + L;
                if (down >= N) down -= N;

                int left = index - 1;
                if (left < 0) left += N;

                int right = index + 1;
                if (right >= N) right -= N;

                int direction[4] = {up, right, down, left}; //can replace 4 with 2*dim

                //CHECK UP DIRECTION
                for (int j = 0; j < 4; j++){
                    if (S[direction[j]] == S[i] && flag[direction[j]] != 1){
                        flag[direction[j]] = 1; //flag
                        stack.push_back(direction[j]);  //stack
                        cluster_size++;                 //increase size of cluster
                    }
                }
            }
            sizes.push_back(cluster_size);
        }
    }
    double average = accumulate(sizes.begin(), sizes.end(), 0.0)/sizes.size();
    return average;
}

double calculate_deviation(double mag[], int size){
    double mean = 0;
    double mean_sq = 0;

    for (int i = 0; i < size; i++){
        mean += mag[i];
        mean_sq += pow(mag[i],2);
    }

    double deviation = mean_sq/size - pow((mean/size),2);
    return deviation;
}

quantity wolff(double T){
    //creates and initialises the spin array of size N = L * L
    int *S;
    S = new int[N];
    init_hot_wolff(S);

    int repeats = 4000;             //the number of steps to be made in the wolf algorithm
    int sample_start = 2000;        //the number of steps where measurements will begin

    int *temp_cluster;              //flagging array
    temp_cluster = new int[N];

    for (int i = 0; i < N ; i++) temp_cluster[i] = 0; //removes all flags

    int index;                      //temporarily stores the current spin index
    int cluster_size=0;             //measures cluster size
    double beta = 1.0/T;            //calculates beta (1/T)
    double prob = 1 - exp(-2*beta); //calculates probability of inclusion

    vector<int> stack;              //
    double magnetisation[repeats];  //array to store all magnetisations
    double sample_magnetisation[repeats-sample_start];      //array to store magnetisasions to be used for analysis

    //sets the magnetisation matrix to zero
    for (int i = 0; i < repeats; i++){
        magnetisation[i] = 0;
    }

    //calculates the initial magnetisation
    for (int i = 0; i < N; i++){
        magnetisation[0] += S[i];
    }

    int cseed;                      //cseed is the random point where a cluser will grow
    double x;                       //x is the probability of cluster growth

    //the wolf algorithm
    for (int i = 1; i < repeats; i++){
        cseed = rand()%N;           //pick random point
        temp_cluster[cseed] = 1;    //flag it
        stack.push_back(cseed);     //add it to the stack
        cluster_size = 1;           //initialise cluster size

        while (stack.size()>0){
            index = stack.back();   //removes last element from the stack
            stack.pop_back();


            //defines all nearest neighbours for clarity
            //uses helical boundary conditions
            int up = index - L;
            if (up < 0) up += N;

            int down = index + L;
            if (down >= N) down -= N;

            int left = index - 1;
            if (left < 0) left += N;

            int right = index + 1;
            if (right >= N) right -= N;

            int direction[4] = {up, right, down, left}; //can replace 4 with 2*dim

            //checks all directions to grow the cluster if the random number x is less then the probability
            for (int j = 0; j < 4; j++){
                x = (double)rand()/RAND_MAX;
                if (x < prob){
                    if (S[direction[j]] == S[cseed] && temp_cluster[direction[j]] != 1){
                        temp_cluster[direction[j]] = 1; //flag
                        stack.push_back(direction[j]);  //stack
                        cluster_size++;                 //increase size of cluster
                    }
                }
            }
        }

        //flips the cluster if its smaller than half the entire size
        if (cluster_size < N/2){
            for (int k = 0; k < N; k++){
                if (temp_cluster[k]){
                    S[k] *= -1;
                    temp_cluster[k] = 0;                //remove flags
                }
            }
            //calculates magnetisation
            magnetisation[i] = (magnetisation[i-1]+cluster_size*S[cseed]*2);
        }
        //flips all but the selected cluster, if the cluster is bigger than half the entire size
        else{
            for (int k = 0; k < N; k++){
                if (!temp_cluster[k]){
                    S[k] *= -1;
                }
                temp_cluster[k] = 0;                   //remove flags
            }
            //calculates magnetisations
            magnetisation[i] = -magnetisation[i-1]+cluster_size*S[cseed]*2;
        }
    }
    //takes the absolute of the magnetisation
    for (int i = 0; i < repeats; i++){
        magnetisation[i] = fabs((magnetisation[i]));
    }

    //pikcs the required range to be used for sampling
    for (int i = 0; i < repeats-sample_start; i++){
        sample_magnetisation[i] = (magnetisation[i+sample_start]/N);
    }

    //defines the autocorrelation array
    double *autocorr;
    autocorr = new double[repeats-sample_start];

    //calculates the autocorrelation function
    autocorr = autocorrelation(sample_magnetisation, repeats-sample_start);

    //finds the 14% mark, i.e. 2 tau_0
    int mark = 0;
    for (int i = 0; i < repeats-sample_start; i++){
        if (autocorr[i] < 0.14){
            mark = i;
            break;
        }
    }

    //temp is used to store the return values from bootstrap (4 values)

    quantity return_value;
    return_value.value = calculate_deviation(sample_magnetisation, repeats-sample_start)*beta*N;    //calculates chi
    return_value.error = bootstrap_error(sample_magnetisation, repeats-sample_start, 512)*beta*N*sqrt(1+mark); //calculates bootstrap error multiplied by sqrt (1+2_tau0)

    //    //output magnetisation and autocorrelation to file if needed
    //    ofstream myfile;
    //    myfile.open("wolf_mag.txt");
    //    for (int i = 0; i < repeats-sample_start; i++){
    //        myfile << i << "," << autocorr[i] <<"," << sample_magnetisation[i]<<endl;
    //    }
    //
    //    if(myfile.is_open()) {
    //        myfile.close();
    //    }

    //delete dynamically allocated arrays
    delete[] S;
    delete[] temp_cluster;
    delete[] autocorr;

    return return_value;
}

void init_hot_wolff(int arr[]){
    for (int i = 0; i < V; i++){
        arr[i] = 2*(rand()%2) - 1;
    }
}

void print_2D_array(int arr[]){
    for (int i = 0; i < V; i++){
        if (arr[i] == 1)
            cout << "1" << " ";
        else if (arr[i] == -1)
            cout << "0" << " ";
        if ((i+1) % L == 0)
            cout << endl;
    }
}
