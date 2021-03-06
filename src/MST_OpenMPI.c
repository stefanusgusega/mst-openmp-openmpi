#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>

#define INF 100005

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

    graph->edges = (struct Edge*)malloc(sizeof(struct Edge)*E);

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
    printf("masuk scatter edges\n");
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype mpiFoo;
    MPI_Type_contiguous( 3, MPI_INT, &mpiFoo );
    MPI_Type_commit( &mpiFoo );
    
    MPI_Scatter(edgeList, *num_edge_part, mpiFoo, edgeListPart, *num_edge_part, mpiFoo, 0, MPI_COMM_WORLD);

    // kalo udah rank terakhir dan gabisa habis dibagi
    // if (rank == size - 1 && num_edge % *num_edge_part != 0) {
    //     *num_edge_part = num_edge%*num_edge_part;
    // }

    // else if (num_edge/2+1 < size && num_edge!=size) {
    //     printf("ada yg mau error rank %d\n",rank);
        
       
    // }
}

void merge(struct Edge* edgeList, const int start, const int mid, const int end) {
    printf("masuk merge\n");
    int left_len = mid - start + 1;
    int right_len = end - mid;
    // printf("nilai right len: %d\n",right_len);
    // printf("nilai end mid: %d %d\n",end, mid);

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

    int i=0,j=0,k = start;
    // printf("nilai k init: %d\n",k);
    while (i < left_len && j < right_len) {
        if (leftside[i].weight <= rightside[j].weight) {
            edgeList[k] = leftside[i++];
        }
        else {
            edgeList[k] = rightside[j++];
        }
        k++;
    }
    // printf("nilai k while pertama: %d\n",k);

    // if there are still any
    while (i < left_len) {
        edgeList[k] = leftside[i++];
        k++;
    }

    while (j < right_len) {
        edgeList[k] = rightside[j++];
        k++;
    }
    // printf("nilai k: %d\n",k);
    for (int it = 0; it < k; it++) {
        printf("merge %d-%d=%d\n",edgeList[it].first,edgeList[it].sec,edgeList[it].weight);
    }

}

void mergeSort(struct Edge* edgeList, int start, int end) {
    int mid = start+(end-start)/2;
    // printf("mergesort %d\n",mid);
    // for (int i = 0; i <= end; i++) {
    //     printf("baru msk merge sort: %d-%d=%d\n",edgeList[i].first,edgeList[i].sec,edgeList[i].weight);
    // }
    // if (start != 2 || end!= 3) {
    //     printf("start end: %d %d\n",start,end);

    // }
    if (start<end){
        mergeSort(edgeList, start, mid);
        // printf("udah masuk merge left\n");
        mergeSort(edgeList, mid+1, end);
        // printf("udah masuk merge right\n");
        merge(edgeList, start, mid, end);

    }

}

void parallelSort(struct Graph* graph) {
    // printf("paralel sort\n");
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
        scatterEdges(graph->edges, edge_list_part, num_edges,&num_edge_part);

    }
    else {
        edge_list_part = graph->edges;
    }

    //sort
    printf("num_edge_part %d: %d\n",rank,num_edge_part);
    for (int k = 0; k < num_edge_part; k++) {
        printf("edges ke %d dari proses %d: %d-%d\n",k,rank,edge_list_part[k].first,edge_list_part[k].sec);
    }
    mergeSort(edge_list_part,0,num_edge_part-1);
    // defining datatype edge
    MPI_Datatype MPI_Edge;
    MPI_Type_contiguous( 3, MPI_INT, &MPI_Edge);
    MPI_Type_commit( &MPI_Edge );
    // merge all sorting result
    if (isParallel) {
        int src,dest,elmt_recv;
        for (int step = 1; step < size; step *= 2) {
            // kalo sedang ada di rank yg kelipatan 2, karena kebagi 2 terus scatternya
            if (rank % (2*step) == 0) {
                src = rank+step;
                // receive how many element
                MPI_Recv(&elmt_recv,1,MPI_INT,src,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                // realloc edge list
                edge_list_part = realloc(edge_list_part,sizeof(struct Edge)*(elmt_recv+num_edge_part));
                // receive the edge list
                MPI_Recv(&edge_list_part[num_edge_part],elmt_recv, MPI_Edge, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
                for (int i = 0; i < num_edges; i++) {
                    printf("setelah reecv:%d-%d=%d\n",edge_list_part[i].first,edge_list_part[i].sec,edge_list_part[i].weight);
                }
                // merge the received list
                merge(edge_list_part, 0, num_edge_part-1, num_edge_part+elmt_recv-1);
                num_edge_part += elmt_recv;
            }
            else {
                dest = rank-step;
                MPI_Send(&num_edge_part,1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&edge_list_part[elmt_recv], elmt_recv, MPI_Edge, dest, 0, MPI_COMM_WORLD);
            }
        }
        if (rank == 0) {
            graph->edges = edge_list_part;
        }
    }
    else {
        graph->edges = edge_list_part;
    }

}

int myComp(const void* a, const void* b)
{
    struct Edge* a1 = (struct Edge*)a;
    struct Edge* b1 = (struct Edge*)b;
    return a1->weight > b1->weight;
}

int comparator(const void* a, const void*b) {
    int i = ((struct Edge*) a)->first;
    int j = ((struct Edge*) b)->first;
    int k = ((struct Edge*) a)->sec;
    int l = ((struct Edge*) b)->sec;
    if (i < j) {
        return -1;
    }
    else {
        if (k < l) {
            return -1;
        }
    }
    return 1;
    
    
}

void kruskal(struct Graph* graph) {
    // masi sequential
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    int num_vertex = graph->num_vertex;
    // printf("num vertex: %d\n",num_vertex);
    struct Edge res[num_vertex];
    // struct Edge* res = (struct Edge*) malloc (sizeof(struct Edge)*num_vertex);
    int res_it=0, edg_it=0;


    //sort
    parallelSort(graph);
    // qsort(graph->edges, graph->num_edge, sizeof(graph->edges[0]),
    //     myComp);
    if (rank==0) {

        struct Subset* subsets = (struct Subset*)malloc(sizeof(struct Subset)*num_vertex);

        // init parent and rank
        for (int ver = 0; ver < num_vertex; ++ver) {
            subsets[ver].parent = ver;
            subsets[ver].rank = 0;
        }
        int minCost = 0;
        while(res_it < num_vertex-1 && edg_it < graph->num_edge) {
            struct Edge next_edge = graph->edges[edg_it++];

            int x = find(subsets, next_edge.first);
            int y = find(subsets, next_edge.sec);

            if (x != y) {
                res[res_it] = next_edge;
                minCost += res[res_it++].weight;
                graph_union(subsets,x,y);
            }

        }
        // gettimeofday(&stop, NULL);

        qsort(&res, res_it, sizeof(struct Edge),comparator);
        printf(
            "%d\n",minCost);
        for ( edg_it= 0; edg_it < res_it; ++edg_it)
        {
            printf("%d-%d\n", res[edg_it].first,
                res[edg_it].sec);
        }
    }
    // printf("Waktu Eksekusi: %lu ms\n", (stop.tv_sec - start.tv_sec) * 1000 + (stop.tv_usec - start.tv_usec) / 1000);
    // return;
}

int main (int argc, char* argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // struct Graph* graph = createGraph(V, E)
    // make the graph
    // struct Graph* graph;
    struct Graph* graph = &(struct Graph ) { .num_edge = 0, .num_vertex = 0, .edges = NULL };
    if (rank == 0) {
        int n, inp;
        int it=0;
        scanf("%d", &n);
        int adj[n][n];
        int V=n; // Number of vertices in graph
        int E=0; // Number of edges in graph
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                scanf("%d", &inp);
                adj[i][j] = inp;
                if(inp != -1 && j >= i) {
                    E++;
                }
            }
        }
        graph = createGraph(V, E);
        for (int i = 0; i < n; i++) {
            for(int j = i; j < n; j++) {
                if (adj[i][j] != -1) {
                    // printf("masuk sini\n");
                    graph->edges[it].first=i;
                    graph->edges[it].sec=j;
                    graph->edges[it].weight=adj[i][j];
                    // printf("kluar\n");
                    it++;
                }
            }
        }
        printf("%d",E);
        
    
        // // add edge 0-1
        // graph->edges[0].first = 0;
        // graph->edges[0].sec= 1;
        // graph->edges[0].weight = 10;
    
        // // add edge 0-2
        // graph->edges[1].first = 0;
        // graph->edges[1].sec = 2;
        // graph->edges[1].weight = 6;
    
        // // add edge 0-3
        // graph->edges[2].first = 0;
        // graph->edges[2].sec= 3;
        // graph->edges[2].weight = 5;
    
        // // add edge 1-3
        // graph->edges[3].first = 1;
        // graph->edges[3].sec= 3;
        // graph->edges[3].weight = 15;
    
        // // add edge 2-3
        // graph->edges[4].first = 2;
        // graph->edges[4].sec= 3;
        // graph->edges[4].weight = 4;

    }
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    kruskal(graph);
    if (rank == 0) {
        gettimeofday(&stop, NULL);
        printf("Waktu Eksekusi: %lu ms\n", (stop.tv_sec - start.tv_sec) * 1000 + (stop.tv_usec - start.tv_usec) / 1000);

    }

    MPI_Finalize();
    return 0;
}