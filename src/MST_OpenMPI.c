#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

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
    int num_edges;
    bool isParallel = size > 1;

    // broadcast jumlah elemen
    if (rank == 0) {
        // ambil jml elemen, bc
        num_edges = graph->num_edge;
        MPI_Bcast(&num_edges,1,MPI_INT,0,MPI_COMM_WORLD);
    }
    else {
        // bc doang
        MPI_Bcast(&num_edges,1,MPI_INT,0,MPI_COMM_WORLD);
    }

    // scatter edgelist
    int num_edge_part = (num_edges+size-1)/size;
    struct Edge* edge_list_part = (struct Edge*)malloc(num_edge_part*sizeof(struct Edge));
    if (isParallel) {
        scatterEdges(graph->edges, edge_list_part, num_edges,num_edge_part);

    }
    else {
        edge_list_part = graph->edges;
    }

    //sort
    mergeSort(&edge_list_part,0,num_edges-1);

    // merge all sorting result
    if (isParallel) {
        int src,dest,elmt_recv;


    }

}

int myComp(const void* a, const void* b)
{
    struct Edge* a1 = (struct Edge*)a;
    struct Edge* b1 = (struct Edge*)b;
    return a1->weight > b1->weight;
}

void kruskal(struct Graph* graph) {
    // masi sequential
    int num_vertex = graph->num_vertex;
    struct Edge res[num_vertex];
    // struct Edge* res = (struct Edge*) malloc (sizeof(struct Edge)*num_vertex);
    int res_it=0, edg_it=0;


    //sort
    // parallelSort(graph);
    qsort(graph->edges, graph->num_edge, sizeof(graph->edges[0]),
        myComp);

    struct Subset* subsets = (struct Subset*)malloc(sizeof(struct Subset)*num_vertex);

    // init parent and rank
    for (int ver = 0; ver < num_vertex; ++ver) {
        subsets[ver].parent = ver;
        subsets[ver].rank = 0;
    }

    while(res_it < num_vertex-1 && edg_it < graph->num_edge) {
        struct Edge next_edge = graph->edges[edg_it++];

        int x = find(subsets, next_edge.first);
        int y = find(subsets, next_edge.sec);

        if (x != y) {
            res[res_it++] = next_edge;
            graph_union(subsets,x,y);
        }

    }
    printf(
        "Following are the edges in the constructed MST\n");
    int minimumCost = 0;
    for ( edg_it= 0; edg_it < res_it; ++edg_it)
    {
        printf("%d -- %d == %d\n", res[edg_it].first,
            res[edg_it].sec, res[edg_it].weight);
        minimumCost += res[edg_it].weight;
    }
    printf("Minimum Cost Spanning tree : %d",minimumCost);
    return;
}

int main (int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int V = 4; // Number of vertices in graph
    int E = 5; // Number of edges in graph
    struct Graph* graph = createGraph(V, E);
 
    // add edge 0-1
    graph->edges[0].first = 0;
    graph->edges[0].sec= 1;
    graph->edges[0].weight = 10;
 
    // add edge 0-2
    graph->edges[1].first = 0;
    graph->edges[1].sec = 2;
    graph->edges[1].weight = 6;
 
    // add edge 0-3
    graph->edges[2].first = 0;
    graph->edges[2].sec= 3;
    graph->edges[2].weight = 5;
 
    // add edge 1-3
    graph->edges[3].first = 1;
    graph->edges[3].sec= 3;
    graph->edges[3].weight = 15;
 
    // add edge 2-3
    graph->edges[4].first = 2;
    graph->edges[4].sec= 3;
    graph->edges[4].weight = 4;

    printf("Hello world!");
    MPI_Finalize();
    return 0;
}