#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

// Define edge
struct Edge {
    int first, sec, weight;
};

// representing graph
struct Graph {
    int num_vertex, num_edge;

    struct Edge* edges;
};

// to show parent of vertex and on top of how many vertex
struct Subset {
    int parent;
    int rank;
};

struct Graph* createGraph(int V, int E) {
    struct Graph* graph = (struct Graph*)malloc(sizeof(struct Graph));
    graph->num_edge = E;
    graph->num_vertex = V;

    graph->edges = (struct Edge*)malloc(sizeof(struct Edge));

    return graph;
};

// to find the root
int find(struct Subset subsets[], int i) {
    if (subsets[i].parent != i) {
        subsets[i].parent = find(subsets, subsets[i].parent);
    }
    return subsets[i].parent;
}

// make a x---y as the member of set
void graph_union(struct Subset subsets[], int x, int y) {
    // seek the root
    int xroot = find(subsets,x);
    int yroot = find(subsets,y);

    // attach the smaller rank tree to the bigger one
    // Union by rank
    if (subsets[xroot].rank < subsets[yroot].rank) {
        subsets[xroot].parent = yroot;
    }
    else if (subsets[xroot].rank > subsets[yroot].rank) {
        subsets[yroot].parent = xroot;
    }
    // if same
    else {
        subsets[yroot].parent = xroot;
        subsets[xroot].rank++;
    }
    
}


void scatterEdges(struct Edge* edgeList, struct Edge* edgeListPart, const int num_edge, int* num_edge_part) {
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype mpiFoo;
    MPI_Type_contiguous( 3, MPI_INT, &mpiFoo );
    MPI_Type_commit( &mpiFoo );
    
    MPI_Scatter(edgeList, *num_edge_part, mpiFoo, edgeListPart, *num_edge_part, mpiFoo, 0, MPI_COMM_WORLD);
}

int comparator(const void* a, const void* b) {
    struct Edge* a1 = (struct Edge*)a;
    struct Edge* b1 = (struct Edge*)b;
    return a1->weight > b1->weight;
}


void merge(struct Edge* edgeList, const int start, const int mid, const int end) {
    int left_len = mid - start + 1;
    int right_len = end - mid;

    // temp arrays
    struct Edge* leftside = (struct Edge*) malloc (left_len*sizeof(struct Edge));
    struct Edge* rightside = (struct Edge*) malloc (right_len*sizeof(struct Edge));

    // pindahin ke temp arrays
    for (int i = 0; i < left_len; i++) {
        leftside[i] = edgeList[start+i];
    }
    for (int j = 0; j < right_len; j++) {
        rightside[j] = edgeList[mid+1+j];
    }

    int i=0,j=0,k = 0;

    while (i < left_len && j < right_len) {
        if (leftside[i].weight <= rightside[j].weight) {
            edgeList[k] = leftside[i++];
        }
        else {
            edgeList[k] = rightside[j++];
        }
        k++;
    }

    // if there are still any
    while (i < left_len) {
        edgeList[k] = leftside[i++];
        k++;
    }

    while (j < right_len) {
        edgeList[k] = rightside[j++];
        k++;
    }

}

void mergeSort(struct Edge* edgeList, int start, int end) {
    int mid = (start+end)/2;
    mergeSort(edgeList, start, mid);
    mergeSort(edgeList, mid+1, end);
    merge(edgeList, start, mid, end);

}

void parallelSort(struct Graph* graph) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    bool isParallel = size > 1;

    // broadcast jumlah elemen
    if (rank == 0) {
        // ambil jml elemen, bc
    }
    else {
        // bc doang
    }

}

void kruskal(struct Graph* graph, struct Graph* mst) {

}

int main (int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    return 0;
}