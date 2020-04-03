#include "header.h"

/*
""Beating critical slowing down requires new algorithms so that at each
step a spin configuration is changed at the scale of a spin cluster"" - K.N.ANAGNOSTOPOULOS
*/

void wolff_function(double Temp, int cycles) {
    /*
    Wolff algorithm:
    1. Choose a random site on the Lattice to begin forming a cluster; the cluster seed.
    2. Repeat: For each new member of the cluster visit its neighbours and if the spin is the same add them to the cluster with probability = 1 - 2 * exp(-2*beta).
    3. When there are no more new members the cluster is complete.
    4. Flip all spin sites of the cluster.
    */

    // useful initialisations
    vector<int> temp_vec; // temporary vector for the nearest neighbours
    double prob = 1 - exp(-2.0/Temp); // probability of adding to cluster
    int cseed; // cluster seed
    double x; // probability variable
    vector<int> new_members; // stack to keep a queue of the new members of the array
    int member; // member index
    int neighbour; // neighbour index

    int *flags; // flag array to keep track of the cluster
    flags = new int[N];
    for (int i = 0; i < N ; i++) flags[i] = 0; // initialise to 0; 1 marks the sites already included in the cluster

    // begin algorithm
    for (int i=0; i<cycles; i++) {
        cseed = rand()%N; // choose a new seed
        flags[cseed] = 1; // flag the seed
        new_members.push_back(cseed); // add the seed

        // create the cluster
        while(new_members.size()>0) {
            // implementing breadth first search so that the cluster grows around the seed
            member = new_members.begin(); // retrieves the first element in the vector
            new_members.erase(members.begin()); // remove first element from new_members
            temp_vec = mapOfNearest[member];
            for (int j=0; j<temp_vec.size(); j++) {
                neighbour = temp_vec[j];
                x = ((double)rand()/RAND_MAX);
                if  (x < prob) && (spins[neighbour] == spins[cseed]) && (flags[neighbour] != 1) {
                    new_members.push_back(neighbour); // add to new_members
                    flags[neighbour] = 1; // flag
                }
            }
        }
        // flip the spins of the cluster and return cluster size
        int cluster_size = flip_cluster(flags,spins[cseed]);

        //

    }




}

int flip_cluster(int flags[], int seed_spin) {
    /*
    flips the spins in the completed cluster, records the size of the cluster and updates the magnetisation of the system
    */
    int cluster_size = 0;
    for (int i = 0; i<N; i++) {
        if (flags[i]) {
            spins[i] = -seed_spin;
            cluster_size++;
        }
    }
    M -= 2*seed_spin*cluster_size;
    return cluster_size;
}
