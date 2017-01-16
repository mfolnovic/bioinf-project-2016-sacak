#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

typedef unsigned char uchar;
typedef unsigned int uint;

// least negative integer
const uint EMPTY = ((uint)1)<<(sizeof(uint)*8-1);

#define scan_rtl_check_S_type(prev, curr) (prev < curr || (prev == curr && current_s_type))

#define scan \
                for (uint i = processing_S_type ? n - 1 : 0; \
                     processing_S_type ? i > 0 : i < n; \
                     processing_S_type ? i -- : i ++)

#define scan_complex(is_processing_S_type, is_processing_LMS, step) \
                for (uint i = is_processing_S_type ? n - 1 : (is_processing_LMS ? n -2 : 0); \
                     is_processing_S_type | is_processing_LMS ? i > 0 : i < n; \
                     is_processing_S_type | is_processing_LMS ? i -= step : i += step)


inline void printArray(void* arr, int n) {
    for (int i = 0; i < n; i++) {
        printf ("%d ", ((int*)arr)[i]);
    }
    printf("\n");
}

inline int shiftValue(int *SA, int start, int direction, bool check_empty) {
    uint count = 0, i;
    // The end item of the left/right-neighbouring bucket can be found by
    // scanning from SA[start] to the left/right, until we see the first item
    // SA[i] that is negative for being reused as a counter.
    for (i = start + direction;
        SA[i] >= 0 || (check_empty && SA[i] == EMPTY);
        i += direction, count += 1);

    int *source, *destination;
    if (direction > 0) {
        source = SA + start;
        destination = source + 1;
    } else {
        source = SA + start - count;
        destination = source - 1;
    }

    // Having found SA[x], we shift-left/right one step all the items in SA[x + 1, c]
    // to SA[x, c-1] and set SA[c] as empty.
    memmove((void*)destination, (const void*)source, sizeof(int) * (count + 1));

    return i;
}

inline int shiftCount(int *SA, int start, int count, int direction) {
    for (int i = 0; i < count; i++) {
        SA[start + i*direction] = SA[start + (i + 1) * direction];
    }
    return count;
}

void initializeBuckets(uchar* T, uint* bkt, uint K, uint n, bool set_to_end) {
    memset(bkt, 0, sizeof(uint) * K);

    // count each character from alphabet
    for (uint i = 0; i < n; i++) {
        bkt[T[i]] += 1;
    }

    uint end_index = -1;
    uint begin_index = 0;
    for (uint i = 0; i < K; i++) {
        int current = bkt[i];

        end_index += current;
        // set index to start/end of the bucket
        bkt[i] = set_to_end ? end_index : begin_index;
        begin_index += current;
    }
}

inline void initializeSA(uint* SA, uint n) {
    memset(SA, 0, sizeof(uint) * n);
}

void fillLMSBuckets0(uchar* T, uint* SA, uint* bkt, uint K, uint n) {
    bool current_s_type = false;

    // initialize buckets to end of each bucket
    initializeBuckets(T, bkt, K, n, /* set_to_end */ true);

    // scan LMS right to left
    for (int i = n-2; i > 0; i--) {
        uchar prev = T[i - 1], curr = T[i];
        bool previous_s_type = scan_rtl_check_S_type(prev, curr);

        // if LMS
        if (!previous_s_type && current_s_type) {
            SA[bkt[curr]] = i;
            bkt[curr] -= 1; // decrease counter
        }

        current_s_type = previous_s_type;
    }

    SA[0] = n - 1;
}

void inducedSort0(uchar* T, uint* SA, uint* bkt, uint K, uint n, bool processing_S_type, bool suffix) {
    // initialize buckets to start/end of each bucket
    initializeBuckets(T, bkt, K, n, /* set_to_end */ processing_S_type);

    if (!processing_S_type) {
        bkt[0] += 1;
    }

    scan {
        // for each scanned non-empty item SA[i]
        if (SA[i] > 0) {
            uint j = SA[i] - 1;
            uchar c, curr, next;
            c = curr = T[j], next = T[j + 1];

            bool is_S_type = curr <= next && bkt[T[j]] < i;
            bool is_L_type = curr >= next;

            // if T[j] is L-type/S-type
            if ((processing_S_type && is_S_type) || (!processing_S_type && is_L_type)) {
                // set SA[bkt[c]] = j
                SA[bkt[curr]] = j;
                // ... and increase/decrease bkt[c] by 1
                bkt[c] += processing_S_type ? -1 : 1;

                if (!suffix && (processing_S_type || (!processing_S_type && i > 0))) SA[i] = 0;
            }
        }
    }
}

/*
  processing_type: 0 - LMS, 1 - S-type, 2 - L-type
*/
void inducedSort1(int* T, int* SA, uint n, int processing_type, bool suffix) {
    int step = 1;

    bool processing_LMS = processing_type == 0;
    bool processing_S_type = processing_type == 1;
    bool processing_L_type = processing_type == 2;
    int mul = processing_LMS || processing_S_type ? -1 : 1;

    if (processing_LMS) {
        // initialize each item of SA[0, N-1] as empty
        for (uint i = 0; i < n; i++) {
            SA[i] = EMPTY;
        }
    }

    bool current_s_type = false;
    scan_complex(processing_S_type, processing_LMS, step) {
        step = 1;
        if (!processing_LMS && SA[i] <= 0) continue;

        uint j;
        int curr;
        bool is_S_type = false, is_L_type = false, is_LMS = false;
        if (processing_LMS) {
            j = i;
            curr = T[j]; int prev = T[j - 1];
            bool previous_s_type = scan_rtl_check_S_type(prev, curr);
            is_LMS = !previous_s_type && current_s_type;
            current_s_type = previous_s_type;
        } else {
            j = SA[i] - 1;
            curr = T[j]; int next = T[j + 1];
            is_S_type = curr < next || (curr == next && curr > (int)i);
            is_L_type = curr >= next;
        }

        if ((processing_LMS && !is_LMS) || (processing_S_type && !is_S_type) || (processing_L_type && !is_L_type)) {
            continue;
        }

        int d = SA[curr];
        // if SA[curr] stores suffix index
        if (d >= 0) {
            // ... then SA[curr] is "borrowed" by the right-neighbouring bucket
            // (of bucket(SA, T, j)). In this case, SA[curr] is storing the smallest
            // item in right-neighbouring bucket, and we need to shift-right one step
            // all the items in the right-neighbouring bucket to their correct locations
            // in SA.
            uint h = shiftValue(SA, curr, -mul, /* check_empty */ !processing_LMS);
            if ((processing_S_type && h > i) || (processing_L_type && h < i)) {
                step = 0;
            }
            d = SA[curr] = EMPTY;
        }

        // if SA[curr] stores empty value...
        if ((uint)d == EMPTY) {
            // ... then suf(T, j) is the first suffix being put into its bucket. In this case
            //, we further check SA[curr +- 1] to see if it is empty or not. If it is...
            if (((processing_LMS || processing_S_type) && (uint)SA[curr - 1] == EMPTY) || (processing_L_type && (uint)curr < n - 1 && (uint)SA[curr + 1] == EMPTY)) {
                // ... we sort suf(T, j) into SA[curr +- 1] by settings SA[curr +- 1] = j and
                // start to reuse SA[curr] as a counter by setting SA[curr] = -1.
                SA[curr + mul] = j;
                SA[curr] = -1;
            } else {
                // Otherwise, SA[curr - 1] may be non-negative for a suffix index or negative
                // for a counter, and suf(T, j) must be the only element of its bucket, we hence
                // simply put suf(T, j) into its bucket by settings SA[curr] = j
                SA[curr] = j;
            }
        } else { // else SA[curr] stores bucket counter
            // In this case, let d = SA[curr] and pos = c - d + 1 / c + d - 1, then SA[pos] is
            // the item that suf(T, j) should be stored into. However, suf(T, j) may be the
            // smallest/largest suffix in its bucket. Therefore, we further check the value
            // of SA[pos] to proceed as follows.
            uint pos = curr - mul * (d - 1);
            if ((uint)SA[pos] == EMPTY && (!processing_L_type || pos <= n-1)) {
                // if SA[pos] is empty, we simply put suf(T, j) into its bucket by setting
                // SA[pos] = j, and increase the counter of its bucket by 1, i.e.
                // SA[curr] = SA[curr] - 1 (notice that SA[curr] is negative for a counter).
                SA[curr] -= 1;
                SA[pos] = j;
            } else {
                // Otherwise, it indicates that SA[pos] is the start item of the
                // right-neighbouring bucket, which must be currently non-negative for a
                // suffix index or negative for a counter. Hence, we need to shift-left/right
                // one step the items in SA[pos + 1, curr - 1]/SA[curr + 1, pos - 1] to
                // SA[pos + 2, curr]/SA[curr, pos - 2], then sort suf(T, j) into
                // its bucket by setting SA[pos -+ 1] = j
                shiftCount(SA, curr, -d, mul);
                SA[pos - mul] = j;
                if ((processing_S_type && T[j] > (int)i) || (processing_L_type && (uint)curr < i)) {
                    step = 0;
                }
            }
        }

        bool is_L_type1 = (j+1 < n-1) && (T[j+1]>T[j+2] || (T[j+1] == T[j+2] && T[j+1] < i));
        if ((processing_L_type && (!suffix || !is_L_type1) && i > 0) || (processing_S_type && !suffix)) {
            SA[step == 0 ? i - mul : i] = EMPTY;
        }
    }

    if (processing_LMS || processing_L_type || !suffix) {
        // scan to shift-right the items in each bucket
        // with its head being reused as a counter
        processing_S_type |= processing_LMS;
        scan {
            int j = SA[i];
            if (j < 0 && (uint)j != EMPTY) { // SA[i] stores bucket counter
                //printf("%d %d %d\n", i, -j, mul, i + -j * mul);
                shiftCount(SA, i, -j, mul);
                SA[i + -j * mul] = EMPTY;
            }
        }
    }

    if (processing_LMS) {
        SA[0] = n - 1;
    }
}

uint problemReduction(uint* SA, uint** T1, uint** SA1, uint n, uint m, int level) {
    uint n1 = 0;

    for (uint i = 0; i < n; i++) {
        if ((!level && SA[i] > 0) || (level && ((int*)SA)[i] > 0)) {
            SA[n1 ++] = SA[i];
        }
    }

    *SA1 = SA;
    *T1 = SA + m - n1;

    return n1;
}

uint computeLexicographicNames(uchar* T, uint* T1, uint* SA, uint n, uint m, uint n1, int level) {
    // init
    for (uint i = n1; i < n; i++) {
        SA[i] = EMPTY;
    }

    // scan SA_1 once from left to right to name each LMS-substring of T by the
    // start position of the substring's bucket in SA_1, resulting in an interim
    // reduced string denoted by Z_1 (where each character points to the start
    // of its bucket in SA_1)
    int previous_lms_length = 0, previous_x = 0;
    int current_name = 0, n_names = 0;
    
    for (uint j = 0; j < n1; j++) {
        int x = SA[j];
        int lms_length = 1;

        // traverse the LMS-substring from its first character T[x] until we
        // see a character T[x + i] less than its preceding T[x + i - 1]. Now,
        // T[x + i - 1] must be L-type.
        uint i;
        for (i = 1; x + i < n && T[x + i] >= T[x + i - 1]; i++);

        // Continue to traverse the remaining characters of the LMS-substring
        // and terminate when we see a character T[x + i] greater than its
        // preceding T[x + i - 1] or T[x + i] is the sentinel. At this point,
        // we know that the start of the succeeding LMS-substring has been
        // traversed and its position was previously recorded when we saw
        // T[x + i] < T[x + i - 1] the last time.
        for (; x + i <= n && T[x + i] <= T[x + i - 1]; i++) {
            if (x + i == n - 1 || T[x + i] > T[x + i - 1]) {
                lms_length = i + 1;
            }
        }

        // now determine if current LMS-substring is different than the last one
        bool is_different = false;
        if (lms_length != previous_lms_length) is_different = true;
        else {
            for (int offset = 0; offset < lms_length && !is_different; offset++) {
                uint current_pos = x + offset;
                uint previous_pos = previous_x + offset;
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
    for (uint i = n - 1, j = m - 1; i >= n1; i--)
        if (SA[i] != EMPTY) SA[j--] = SA[i];

    // Scan Z_1 once from right to left to replace each S-type character in Z1
    // by the end position of its bucket in SA_1, resulting in the new string T1.
    bool current_s_type = false;
    for (int i = n1-1; i > 0; i--) {
        bool previous_s_type = scan_rtl_check_S_type(T1[i - 1], T1[i]);

        // if S-type
        if (previous_s_type) {
            T1[i - 1] += SA[T1[i - 1]] - 1;
        }

        current_s_type = previous_s_type;
    }

    // return number of unique LMS substrings (or number of names...)
    return n_names;
}

// finds LMS positions in T and stores it in T1 (in text order)
void getSALMS(uint *SA, uchar *T, uint *T1, uint n, uint n1, int level) {
    uint j = n1 - 1;
    T1[j--] = n - 1;
    bool current_s_type = false;

    for (int i = n-2; i > 0; i--) {
        bool previous_s_type = scan_rtl_check_S_type(T[i - 1], T[i]);

        // if LMS
        if (!previous_s_type && current_s_type) {
            T1[j--] = i;
        }
        current_s_type = previous_s_type;
    }

    for (int i = 0; i < n1; i++) {
        SA[i] = T1[SA[i]];
    }

    // re-initialize SA[n1..n-1]
    for (int i = n1; i < n; i++) {
        SA[i] = level ? EMPTY : 0;
    }
}

void putSuffix0(uint *SA, uchar *T, uint *bkt, uint n, uint K, int n1) {
    // find the ends of each bucket
    initializeBuckets(T, bkt, K, n, true);

    // put the suffixes into their buckets
    for (int i = n1 - 1; i > 0; i--) {
        uint j = SA[i];
        SA[i] = 0;
        SA[bkt[T[j]]--] = j;
    }
    SA[0] = n - 1; // set the sentinel
}

void putSuffix1(int *SA, int *T, int n1) {
    int pos, curr, prev = -1;

    for(int i = n1 - 1; i > 0; i--) {
        uint j = SA[i];
        SA[i] = EMPTY;
        curr = T[j];

        if (curr != prev) {
          prev = curr; pos = curr;
        }

        SA[pos--] = j;
    }
}

void sacak(uchar* T, uint* SA, uint K, uint n, uint m, int level) {
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
        fillLMSBuckets0(T, SA, bkt, K, n);

        // compute into bkt[0, K-1] the start position of each bucket in SA. Scan
        // SA once from left to right to induces sort the L-type suffixes of T
        // into their buckets in SA, from the start to the end in each bucket.
        inducedSort0(T, SA, bkt, K, n, /* processing_S_type */ false, /* suffix */ false);

        // compute into bkt[0, K-1] the end position of each bucket in SA. Scan
        // SA once from right to left to induces sort the S-type suffixes of T
        // into their buckets in SA, from the end to the start in each bucket.
        inducedSort0(T, SA, bkt, K, n, /* processing_S_type */ true, /* suffix */ false);
    } else {
        // induced sort all the LMS-substrings of T, reusing the start or
        // the end of each bucket as the bucket's counter
        inducedSort1((int*)T, (int*)SA, n, /* processing_type */ 0, /* suffix */ false);
        inducedSort1((int*)T, (int*)SA, n, /* processing_type */ 2, /* suffix */ false);
        inducedSort1((int*)T, (int*)SA, n, /* processing_type */ 1, /* suffix */ false);
    }

    // SA is reused for storing T1 and SA1
    uint* T1 = NULL;
    uint* SA1 = NULL;
    uint n1 = problemReduction(SA, &T1, &SA1, n, m, level);
    uint K1 = computeLexicographicNames(T, T1, SA, n, m, n1, level);

    // stage 3: sort recursively
    if (K1 == n1) {
        // Directly compute SA(T1) from T1
        for (uint i = 0; i < n1; i++) SA1[T1[i]] = i;
    } else {
        sacak((uchar*) T1, SA1, K, n1, m - n1, level + 1);
    }
    
    // stage 4: induced sort SA(T) from SA(T1).
    getSALMS(SA, T, T1, n, n1, level);

    if (level == 0) {      
        // Induced sort SA(T) from SA(T1) using bkt for bucket counters
        putSuffix0(SA, T, bkt, n, K, n1);
        inducedSort0(T, SA, bkt, K, n, /* processing_S_type */ false, /* suffix */ true);
        inducedSort0(T, SA, bkt, K, n, /* processing_S_type */ true, /* suffix */ true);
        
        // Free the space allocated for bkt
        free(bkt);
    }
    else {
        // Induced sort SA(T) from SA(T1), reusing the start or the end
        // of each bucket as the bucket's counter
        putSuffix1((int *)SA, (int *)T, n1);
        inducedSort1((int*)T, (int*)SA, n, /* processing_type */ 2, /* suffix */ true);
        inducedSort1((int*)T, (int*)SA, n, /* processing_type */ 1, /* suffix */ true);
    }
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

    clock_t start = clock();
    sacak(T, SA, 256, n, n, 0);
    printArray(SA, n);
    
    double duration = (double)(clock() - start) / CLOCKS_PER_SEC;
    fprintf(stderr, "\nSize: %u bytes, Time: %5.3f seconds\n", n - 1, duration);

    // free allocated memory
    free(T);
    free(SA);

    return 0;
}
