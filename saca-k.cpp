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

void sacak(char* T, int* SA, int K, int n, int level) {
    int *bkt = NULL;

    if (level == 0) {
        bkt = new int[K];

        // initialize buckets to end of each bucket
        initializeBuckets(T, bkt, K, n);

        printf("%d\n", K);
        for (int i = 0; i < K; i++) {
            printf("%d ", *(bkt + i));
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
