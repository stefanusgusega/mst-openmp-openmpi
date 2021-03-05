#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
struct Compare { int val; int index; };    
#pragma omp declare reduction(minimum : struct Compare : omp_out = omp_in.val < omp_out.val ? omp_in : omp_out)
int main()
{
    struct Compare min; 
    int theArray[5];
    theArray[0] = 1;
    theArray[1] = 10;
    theArray[2] = -1;
    theArray[3] = 21;
    theArray[4] = 4;
    min.val = theArray[0]; 
    min.index = 0;
    #pragma omp parallel for reduction(minimum:min)
    for(int i = 1; i<5; i++) {
        if(theArray[i]<min.val) { 
            min.val = theArray[i];
            min.index = i;
        }
    }

    printf("%d %d\n", min.val, min.index);
}

// n = 8
// a = 3, b = 5

// 3+(3+5)+(3+5*2)+...+(3+5*7)

// 0 1 2 3 4 5 6 7 8 9

// 0 3 6 9
