#define LAB4_EXTEND

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab4_IO.h"
#include <mpi.h>
#define EPSILON 0.00001
#define DAMPING_FACTOR 0.85

#define THRESHOLD 0.0001

// void load_data()
// {
// 	int size;
// 	FILE *file_pointer;
// 	int *edges_count_array;
// 	if ((file_pointer = fopen("data_input","r")) == NULL){
// 		exit(1);
// 	}
// 	fscanf(file_pointer,"%d\n",number_of_nodes);

// 	edges_count_array = calloc(number_of_nodes*sizeof(int));
	
// 	int src, dest;
// 	while (!feof(file_pointer)){
// 		f
// 	}


// }


int main (int argc, char* argv[]){
    struct node *nodehead;
    int nodecount;
    int *num_in_links, *num_out_links;
    double *r, *r_pre;
    int i, j;
    double damp_const;
    int iterationcount = 0;
    int collected_nodecount;
    double *collected_r;
    double cst_addapted_threshold;
    double error;
    FILE *fp;
    MPI_Init(&argc,&argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD,&size);

   
	    	// Load the data and simple verification

	     if (get_node_stat(&nodecount, &num_in_links, &num_out_links)) return 254;

	    // Adjust the threshold according to the problem size
	    cst_addapted_threshold = THRESHOLD;
	    
	    // Calculate the result
	    if (node_init(&nodehead, num_in_links, num_out_links, 0, nodecount)) return 254;
	    
	    r = malloc(nodecount * sizeof(double));
	    r_pre = malloc(nodecount * sizeof(double));
	    for ( i = 0; i < nodecount; ++i)
	        r[i] = 1.0 / nodecount;
	    damp_const = (1.0 - DAMPING_FACTOR) / nodecount;


  	int partition = nodecount/size;
  	int remainder = nodecount % size;
    int start = rank*partition;
    int end = rank*partition+partition;
    int *disps = malloc(size*sizeof(int));
    int *partitions = malloc(size*sizeof(int));
    int d;
    for (d=0;d<size;d++){
    	disps[d] = d*partition;
    }

    for (d=0;d<size;d++){
        partitions[d] = partition;
    }
    
    // CORE CALCULATION
    do{
        ++iterationcount;
        vec_cp(r, r_pre, nodecount);
        if (end+partition > nodecount){
        	end = end + remainder;
        }
        for ( i = start; i < end; ++i){
            r[i] = 0;
            for ( j = 0; j < nodehead[i].num_in_links; ++j)
                r[i] += r_pre[nodehead[i].inlinks[j]] / num_out_links[nodehead[i].inlinks[j]];
            r[i] *= DAMPING_FACTOR;
            r[i] += damp_const;
       }
	MPI_Allgatherv(MPI_IN_PLACE,0,MPI_DOUBLE,r,partitions,disps,MPI_DOUBLE,MPI_COMM_WORLD);
    }while(rel_error(r, r_pre, nodecount) >= EPSILON);
    //printf("Program converges at %d th iteration.\n", iterationcount);

    MPI_Finalize();
    Lab4_saveoutput(r,nodecount,0);
    // post processing
    node_destroy(nodehead, nodecount);
    free(num_in_links); free(num_out_links);

}





