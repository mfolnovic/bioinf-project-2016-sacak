#include <cstdio>
#include <cstdlib>
#include <cstring>

void initializeBuckets(char* T, int* bkt, int K, int n) {
    memset(bkt, 0, sizeof(int) * K);

    for (int i = 0; i < n; i++) {
        bkt[T[i]] += 1;
    }

    int end_index = -1;
    int begin_index = 0;
    for (int i = 0; i < K; i++) {
        int current = bkt[i];

        end_index += current;
        bkt[i] = end_index;
        begin_index += current;
    }
}

void initializeSA(int* SA, int n) {
    memset(SA, 0, sizeof(int) * n);
}

void induceSortLMS(char* T, int* SA, int* bkt, int K, int n) {
    bool current_s_type = false;

    for (int i = n-2; i > 0; i--) {
        bool previous_s_type = T[i - 1] < T[i] || (T[i - 1] == T[i] && current_s_type);

        // if LMS
        if (!previous_s_type && current_s_type) {
            SA[bkt[T[i]]] = i;
            bkt[T[i]] -= 1;
        }

        current_s_type = previous_s_type;
    }

    SA[0] = n - 1;
}

void sacak(char* T, int* SA, int K, int n, int level) {
    int *bkt = NULL;

    if (level == 0) {
        bkt = new int[K];

        // initialize buckets to end of each bucket
        initializeBuckets(T, bkt, K, n);

        // initialize each item of SA[0, N-1] as empty
        initializeSA(SA, n);

        // compute into bkt[0, K-1] the end position of each bucket in SA. Scan
        // SA[0, N-1] once from right to left to put all the sorted LMS-suffixes
        // of T into their buckets in SA, from the end to the start in each
        // bucket.
        induceSortLMS(T, SA, bkt, K, n);

        printf("%d\n", K);
        for (int i = 0; i < K; i++) {
            printf("%d ", *(bkt + i));
        }
        printf("\n\n");

        for (int i = 0; i < n; i++) {
            printf("%d ", SA[i]);
        }
        printf("\n");
    } else {

    }

    free(bkt);
}

int main(int argc, char** argv) {
    // open file
    FILE *in = fopen(argv[1], "rb");
   
    // determine file length
    fseek(in, 0, SEEK_END);
    int n = ftell(in) + 1; // +1 because of sentinel
    fseek(in, 0, SEEK_SET);

    char* T = new char[n]; // allocate input string array
    int* SA = new int[n]; // allocate suffix array
    fread(T, sizeof(char), n - 1, in); // read input string
    T[n - 1] = 0; // sentinel

    sacak(T, SA, 256, n, 0);

    // free allocated memory
    free(T);
    free(SA);

    return 0;
}
