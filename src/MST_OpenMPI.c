#include <stdio.h>
#include <mpi.h>
#include <string.h>

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