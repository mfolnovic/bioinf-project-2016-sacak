#include <cstdio>
#include <cstdlib>
#include <cstring>

typedef unsigned char uchar;
typedef unsigned int uint;

// least negative integer
const uint EMPTY = ((uint)1)<<(sizeof(uint)*8-1);

void initializeBuckets(uchar* T, uint* bkt, uint K, uint n, bool set_to_end) {
    memset(bkt, 0, sizeof(uint) * K);

    for (int i = 0; i < n; i++) {
        bkt[T[i]] += 1;
    }

    uint end_index = -1;
    uint begin_index = 0;
    for (int i = 0; i < K; i++) {
        int current = bkt[i];

        end_index += current;
        bkt[i] = set_to_end ? end_index : begin_index;
        begin_index += current;
    }
}

void initializeSA(uint* SA, uint n) {
    memset(SA, 0, sizeof(uint) * n);
}

void fillLMSBuckets(uchar* T, uint* SA, uint* bkt, uint K, uint n) {
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

void inducedSort(uchar* T, uint* SA, uint* bkt, uint K, uint n, bool processing_S_type) {
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
            uint j = SA[i] - 1;
            uchar c = T[j];

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

uint problemReduction(uint* SA, uint** T1, uint** SA1, uint n) {
  uint n1 = 0;

  for (int i = 0; i < n; i++) {
    if (SA[i] > 0) {
        SA[n1 ++] = SA[i];
    }
  }

  *SA1 = SA;
  *T1 = SA + n - n1;

  return n1;
}

uint computeLexicographicNames(uchar* T, uint* T1, uint* SA, uint n, uint n1) {
    // init
    for (int i = n1; i < n; i++) {
        SA[i] = EMPTY;
    }

    // scan SA_1 once from left to right to name each LMS-substring of T by the
    // start position of the substring's bucket in SA_1, resulting in an interim
    // reduced string denoted by Z_1 (where each character points to the start
    // of its bucket in SA_1)
    int previous_lms_length = 0, previous_x = 0;
    int current_name = 0, n_names = 0;
    for (int j = 0; j < n1; j++) {
        int x = SA[j];
        int lms_length = 1;

        // traverse the LMS-substring from its first character T[x] until we
        // see a character T[x + i] less than its preceding T[x + i - 1]. Now,
        // T[x + i - 1] must be L-type.
        int i;
        for (i = 1; x + i < n && T[x + i] >= T[x + i - 1]; i++);

        // Continue to traverse the remaining characters of the LMS-substring
        // and terminate when we see a character T[x + i] greater than its
        // preceding T[x + i - 1] or T[x + i] is the sentinel. At this point,
        // we know that the start of the succeeding LMS-substring has been
        // traversed and its position was previously recorded when we saw
        // T[x + i] < T[x + i - 1] the last time.
        for (; x + i < n && T[x + i] <= T[x + i - 1]; i++) {
            if (x + i == n - 1 || T[x + i] > T[x + i - 1]) {
                lms_length = i + 1;
            }
        }

        // now determine if current LMS-substring is different than the last one
        bool is_different = false;
        if (lms_length != previous_lms_length) is_different = true;
        else {
            for (int offset = 0; offset < lms_length && !is_different; offset++) {
                int current_pos = x + offset;
                int previous_pos = previous_x + offset;
                is_different = current_pos == n - 1 || previous_pos == n - 1 ||
                    T[current_pos] != T[previous_pos];
            }
        }

        if (is_different) {
            // it's different so create new name...
            current_name = j;
            n_names += 1;
            SA[current_name] = 1;

            // it's different so we'll compare next LMS-substring with this one
            previous_x = x;
            previous_lms_length = lms_length;
        } else {
            // it's same, so reuse the name
            SA[current_name] += 1;
        }

        SA[n1 + x/2] = current_name;
    }

    // compact...
    for (int i = n - 1, j = n - 1; i >= n1; i--)
        if (SA[i] != EMPTY) SA[j--] = SA[i];

    // Scan Z_1 once from right to left to replace each S-type character in Z1
    // by the end position of its bucket in SA_1, resulting in the new string T1.
    bool current_s_type = false;
    for (int i = n1-1; i > 0; i--) {
        bool previous_s_type = T[i - 1] < T[i] || (T[i - 1] == T[i] && current_s_type);

        // if S-type
        if (previous_s_type) {
            T1[i - 1] += SA[T1[i - 1]] - 1;
        }

        current_s_type = previous_s_type;
    }

    // return number of unique LMS substrings (or number of names...)
    return n_names;
}

void sacak(uchar* T, uint* SA, uint K, uint n, int level) {
    uint *bkt = NULL;

    // Stage 1: induces sort the LMS-substrings of T
    if (level == 0) {
        bkt = new uint[K];

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

    // SA is reused for storing T1 and SA1
    uint* T1 = NULL;
    uint* SA1 = NULL;
    uint n1 = problemReduction(SA, &T1, &SA1, n);
    uint K1 = computeLexicographicNames(T, T1, SA, n, n1);

    free(bkt);
}

int main(int argc, char** argv) {
    // open file
    FILE *in = fopen(argv[1], "rb");
   
    // determine file length
    fseek(in, 0, SEEK_END);
    int n = ftell(in) + 1; // +1 because of sentinel
    fseek(in, 0, SEEK_SET);

    uchar* T = new uchar[n]; // allocate input string array
    uint* SA = new uint[n]; // allocate suffix array
    fread(T, sizeof(uchar), n - 1, in); // read input string
    T[n - 1] = 0; // sentinel

    sacak(T, SA, 256, n, 0);

    // free allocated memory
    free(T);
    free(SA);

    return 0;
}
