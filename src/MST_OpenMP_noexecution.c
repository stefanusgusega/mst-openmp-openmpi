#include <omp.h>
#include <stdio.h>
#include <sys/time.h>

#define N 3005
#define INF 100005

struct user_definedMin { 
    int val; 
    int index; 
}; 
#pragma omp declare reduction(minimum : struct user_definedMin : omp_out = omp_in.val < omp_out.val ? omp_in : omp_out)

int n, adj[N][N], selected[N], min_edge[N][2], result[N][2];

int tmin(int a, int b) {
    return (a < b) ? a : b;
}

int tmax(int a, int b) {
    return (a > b) ? a : b;
}


int main() 
{
    int threads = omp_get_max_threads();

    scanf("%d", &n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            scanf("%d", &adj[i][j]);
            if(adj[i][j] == -1) adj[i][j] = INF;
        }

        selected[i] = 0;
        min_edge[i][0] = INF;
        min_edge[i][1] = -1;
    }

    int total_weight = 0;
    min_edge[0][0] = 0;

    int cur = 0;

    struct timeval stop, start;
    gettimeofday(&start, NULL);
    for(int i = 0; i < n; i++) {
        int t = -1;

        int index = -1;
        int fmin = INF;
        int j;
        #pragma omp parallel num_threads(threads)
        {
            int index_local = index;
            int fmin_local = fmin;
            #pragma omp for nowait 
            for(j = 0; j < n; j++) {
                if(selected[j] == 0 && min_edge[j][0] < fmin_local) {
                    fmin_local = min_edge[j][0];
                    index_local = j;
                }
            }
            #pragma omp critical 
            {
                if(fmin_local < fmin) {
                    fmin = fmin_local;
                    index = index_local;
                }
            }
        }
        t = index;
        selected[t] = 1;
        total_weight += min_edge[t][0];
        if(min_edge[t][1] != -1) {
            result[cur][0] = tmin(t, min_edge[t][1]);
            result[cur][1] = tmax(t, min_edge[t][1]);
            cur++;
        }
        int to;
        #pragma omp parallel for num_threads(threads) private(to) shared(n)
        for(to = 0; to < n; to++) {
            if(adj[t][to] < min_edge[to][0]) {
                min_edge[to][0] = adj[t][to];
                min_edge[to][1] = t;
            }
        }
    }
    gettimeofday(&stop, NULL);

    for(int i = cur - 1; i >= 1; i--) {
        for(int j = 0; j < i; j++) {
            if(result[j][0] > result[j + 1][0] || (result[j][0] == result[j + 1][0] && result[j][1] > result[j + 1][1])) {
                int temp1 = result[j + 1][0];
                result[j + 1][0] = result[j][0];
                result[j][0] = temp1;

                int temp2 = result[j + 1][1];
                result[j + 1][1] = result[j][1];
                result[j][1] = temp2;
            }
        }
    }

    printf("%d\n", total_weight);

    for(int i = 0; i < cur; i++) {
        printf("%d-%d\n", result[i][0], result[i][1]);
    }



    return 0;
}