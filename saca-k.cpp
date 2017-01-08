#include <cstdio>
#include <cstdlib>
#include <cstring>

void initializeBuckets(char* T, int* bkt, int K, int n, bool set_to_end) {
    memset(bkt, 0, sizeof(int) * K);

    for (int i = 0; i < n; i++) {
        bkt[T[i]] += 1;
    }

    int end_index = -1;
    int begin_index = 0;
    for (int i = 0; i < K; i++) {
        int current = bkt[i];

        end_index += current;
        bkt[i] = set_to_end ? end_index : begin_index;
        begin_index += current;
    }
}

void initializeSA(int* SA, int n) {
    memset(SA, 0, sizeof(int) * n);
}

void fillLMSBuckets(char* T, int* SA, int* bkt, int K, int n) {
    bool current_s_type = false;

    // initialize buckets to end of each bucket
    initializeBuckets(T, bkt, K, n, /* set_to_end */ true);

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

void inducedSort(char* T, int* SA, int* bkt, int K, int n, bool processing_S_type) {
    // initialize buckets to start/end of each bucket
    initializeBuckets(T, bkt, K, n, /* set_to_end */ processing_S_type);

    if (!processing_S_type) {
        bkt[0] += 1;
    }

    for (int i = processing_S_type ? n - 1 : 0;
         processing_S_type ? i > 0 : i < n;
         processing_S_type ? i-- : i++) {
        // for each scanned non-empty item SA[i]
        if (SA[i] > 0) {
            int j = SA[i] - 1;
            char c = T[j];

            bool is_S_type = T[j] <= T[j + 1] && bkt[T[j]] < i;
            bool is_L_type = T[j] >= T[j + 1];

            // if T[j] is L-type/S-type
            if ((processing_S_type && is_S_type) || (!processing_S_type && is_L_type)) {
                // set SA[bkt[c]] = j
                SA[bkt[c]] = j;
                // ... and increase/decrease bkt[c] by 1
                bkt[c] += processing_S_type ? -1 : 1;

                if (i > 0) SA[i] = 0;
            }
        }
    }
}

void sacak(char* T, int* SA, int K, int n, int level) {
    int *bkt = NULL;

    if (level == 0) {
        bkt = new int[K];

        // initialize each item of SA[0, N-1] as empty
        initializeSA(SA, n);

        // compute into bkt[0, K-1] the end position of each bucket in SA. Scan
        // SA[0, N-1] once from right to left to put all the sorted LMS-suffixes
        // of T into their buckets in SA, from the end to the start in each
        // bucket.
        fillLMSBuckets(T, SA, bkt, K, n);

        // compute into bkt[0, K-1] the start position of each bucket in SA. Scan
        // SA once from left to right to induces sort the L-type suffixes of T
        // into their buckets in SA, from the start to the end in each bucket.
        inducedSort(T, SA, bkt, K, n, /* processing_S_type */ false);

        // compute into bkt[0, K-1] the end position of each bucket in SA. Scan
        // SA once from right to left to induces sort the S-type suffixes of T
        // into their buckets in SA, from the end to the start in each bucket.
        inducedSort(T, SA, bkt, K, n, /* processing_S_type */ true); 
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
