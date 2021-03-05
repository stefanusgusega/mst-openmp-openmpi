#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

struct Foo {
    int a,b,c;
};

int main(int argc, char** argv) {
    struct Foo* foo_inst = (struct Foo*) malloc (sizeof(struct Foo)*8);
    for (int i = 0; i < 8; i++) {
        foo_inst[i].a = i;
        foo_inst[i].b = i+2;
        foo_inst[i].c = i*2;
    }
    for (int j = 0; j < 4; j++) {
        printf("%d %d %d\n",foo_inst[j].a, foo_inst[j].b, foo_inst[j].c);
    }
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    MPI_Datatype mpiFoo;
    MPI_Type_contiguous( 3, MPI_INT, &mpiFoo );
    MPI_Type_commit( &mpiFoo );

    struct Foo* container = (struct Foo*) malloc (sizeof(struct Foo)*8);
    MPI_Scatter(foo_inst, 2, mpiFoo, container, 2, mpiFoo, 0, MPI_COMM_WORLD);
    
    printf("isinya %d %d %d %d\n", rank, container[0].a, container[0].b, container[0].c);
    // if (rank != 0) {
    //     MPI_Recv(&foo_inst, 1, mpiFoo, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // }

    if (rank != 0) {
        
    }
    MPI_Finalize();

    return 0;
}
