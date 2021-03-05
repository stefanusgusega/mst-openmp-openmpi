#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

void Trap(int a, int b, int n, int *global_result_p);

int tmin(int a, int b) {
    return (a > b) ? b : a;
}

int main(int argc, char* argv[]) {

    int global_sum = 0;
    int a = 3;
    int b = 5;
    int thread_count = strtol(argv[1], NULL, 10);
    int n = 11;

    #pragma omp parallel num_threads(thread_count)
    {
        Trap(a, b, n, &global_sum);
    }

    int glob = 0;
    #pragma omp parallel num_threads(thread_count)
    {
        Trap(5, 7, 25, &glob);
    }
    

    printf("Result: %d\n", global_sum);
    return 0;
}

// n = 8
// a = 3, b = 5

// 3+(3+5)+(3+5*2)+...+(3+5*7)

// 0 1 2 3 4 5 6 7 8 9

// 0 3 6 9

void Trap(int a, int b, int n, int *global_result_p) {
    int my_rank = omp_get_thread_num();
    int thread_count = omp_get_num_threads();

    int local_n = n / thread_count + !(n % thread_count == 0);

    
    int local_i = local_n * my_rank;
    if(local_i >= n) return;
    int local_max_i = tmin(n - 1, local_i + local_n - 1);

    printf("from %d: %d -> %d\n", my_rank, local_i, local_max_i);

    int local_sum = 0;
    for(int i = local_i; i <= local_max_i; i++) {
        local_sum += a + b * i;
    }

    #pragma omp critical
    *global_result_p += local_sum;
    
}