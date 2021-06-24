#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include <time.h>
#include <m4ri/m4ri.h>

#include "sbox32.h"

typedef struct _crwsm {
    uint8_t row_weight;
    int nrows;
    int ncols;
    uint8_t **rows;
} ConstantRowWeightSparseMatrix;


typedef struct _list {
    int *values;
    int n;
} List;

typedef uint16_t Label;


#define MAX_N 42
#define MAX_L 11
#define MAX_M 15
#define MAX_VPUB_SIGNATURES 1 << MAX_L


#define SBOXSIZE (1 << (SBOX_LOG2_SIZE))

#if defined USE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS && USE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS == 1
    #define MAX_LABEL (1 << (MAX_M + 1))
    #define SBOX_VAL(X) (SBOX[(X % (SBOXSIZE - 1))]);
#else
    #define MAX_LABEL (1 << (SBOX_LOG2_SIZE + 1))
    #define SBOX_VAL(X) (SBOX[(X)])
#endif

typedef uint32_t Sign;

typedef struct _perm_extractor {
    mzd_t* A;
    mzd_t* Vpub;
    mzd_t* Alow;
    mzd_t* A2_inv_times_A1;
    mzd_t* A1;
    int repetitions[MAX_LABEL];
    int n_a1;
    int* A1_indexes;
    int n_a2;
    int* A2_indexes;
    Label* Alabels;
    int n_unique_labels;
    Label* unique_labels;
    List labels_to_Aindexes[MAX_LABEL];
    List labels_to_A1indexes[MAX_LABEL];
    int Vpub_signatures_count[MAX_VPUB_SIGNATURES];
    int Vpub_signatures[MAX_N];
} PermutationExtractor;


typedef struct _poss {
    int* nposs;
    int** poss;
    int nlevels;
    int* support;
    int level;
    ConstantRowWeightSparseMatrix* target_matrix;
    ConstantRowWeightSparseMatrix* pattern_matrix;
    Sign signatures[MAX_M][MAX_M];
    int limits_of_removed_per_level[MAX_M][MAX_M];
} PossibleChildrenPerLevel;


typedef struct _rows_per_distance {
    int* n_for_distance;
    int** data;
    int n;
} RowsPerDistance;


void poss_per_level_add_row_to_support(PossibleChildrenPerLevel* poss_per_level, int row);
void extract_permutation_from_matching(PermutationExtractor* e, ConstantRowWeightSparseMatrix* M, int* rowset, int n_rowset);


Sign get_signature_for_rows(ConstantRowWeightSparseMatrix *matrix, int arr[], int n) {

    uint16_t sigs[matrix->ncols];
    memset(sigs, 0, matrix->ncols*sizeof(*sigs));

    for (int ri = 0; ri < n; ri++) {
        int row = arr[ri];
        for (int ci = 0; ci < matrix->row_weight; ci++) {
            int col = matrix->rows[row][ci];
            sigs[col] += (1 << ri);
        }
    }

    Sign key = 0;

    for (int j = 0; j < matrix->ncols; j++) {
        uint16_t val = sigs[j];
        key += SBOX_VAL(val);
    }
    return key;
}


Sign get_signature_cummulative_base(ConstantRowWeightSparseMatrix *matrix, int arr[], int n, uint16_t sigs[]) {

    memset(sigs, 0, matrix->ncols*sizeof(*sigs));

    for (int ri = 0; ri < n; ri++) {
        int row = arr[ri];
        for (int ci = 0; ci < matrix->row_weight; ci++) {
            int col = matrix->rows[row][ci];
            sigs[col] += (1 << ri);
        }
    }

    Sign key = 0;

    for (int j = 0; j < matrix->ncols; j++) {
        uint16_t val = sigs[j];
        key += SBOX_VAL(val);
    }
    return key;
}


Sign get_signature_cummulative(ConstantRowWeightSparseMatrix *matrix, int new_row, int level, uint16_t sigs[]) {

    for (int ci = 0; ci < matrix->row_weight; ci++) {
        int col = matrix->rows[new_row][ci];
        sigs[col] += (1 << level);
    }

    Sign key = 0;

    for (int j = 0; j < matrix->ncols; j++) {
        // __builtin_prefetch (&(SBOX_VAL(sigs[j])), 0, 1);
        uint16_t val = sigs[j];
        // printf("sigs[%d] = %d\n", j, val);
        key += SBOX_VAL(val);
    }
    return key;
}


ConstantRowWeightSparseMatrix* read_sparse_matrix_input(char filepath[], int max_nrows) {
    ConstantRowWeightSparseMatrix* matrix = malloc(sizeof(*matrix));

    FILE* file = fopen(filepath, "r");
    if (!file) {
        fprintf(stderr, "Could not open file %s.\n", filepath);
        exit(1);
    }
    int total_rows;
    fscanf(file, "%d %d", &total_rows, &matrix->ncols);

    matrix->nrows = max_nrows;
    if (max_nrows == -1) {
        matrix->nrows = total_rows;
    }

    matrix->rows = malloc(matrix->nrows*sizeof(*(matrix->rows)));

    int matrix_weight;
    fscanf(file, "%d", &matrix_weight);

    matrix->row_weight = matrix_weight/total_rows;

    for (int i = 0; i < matrix->nrows; i++) {
        matrix->rows[i] = malloc(matrix->row_weight*sizeof(**(matrix->rows)));

        int this_row_weight = 0;
        for (int j = 0; j < matrix->row_weight; j++) {
            int row;
            int col;
            fscanf(file, "%d %d", &row, &col);
            matrix->rows[i][this_row_weight++] = col;
        }
    }

    fclose(file);
    return matrix;
}

void read_secret_support_file(int secret_support[], int nsecret_support, char filepath[]) {
    FILE *file;
    file = fopen(filepath, "r");
    if (!file) {
        fprintf(stderr, "Could not open file %s.\n", filepath);
        exit(1);
    }

    int n;
    fscanf(file, "%d", &n);
    // printf("read size %d\n", n);

    for (int i = 0; (i < nsecret_support); i++) {
        fscanf(file, "%d", &n);
        secret_support[i] = n;
        // printf("read %d\n", n);
    }

    fclose(file);
}


void test_read_sparse_matrix_input(char filepath[]) {
    ConstantRowWeightSparseMatrix* matrix = read_sparse_matrix_input(filepath, 10);

    printf("nrows, ncols = %d, %d\n", matrix->nrows, matrix->ncols);
    printf("row_weight = %d\n", matrix->row_weight);
    for (int i = 0; i < matrix->nrows; i++) {
        printf("row[%d]: ", i);
        for (int j = 0; j < matrix->row_weight; j++) {
            printf("%3d ", matrix->rows[i][j]);
        }
        printf("\n");
    }

    int list_of_rows[] = {0, 1, 2, 3, 4};
    int nlist_of_rows = sizeof(list_of_rows)/sizeof(list_of_rows[0]);
    printf("nlist_of_rows: %d, matrix->nrows: %d\n", nlist_of_rows, matrix->nrows);
    assert(nlist_of_rows < matrix->nrows);

    get_signature_for_rows(matrix, list_of_rows, nlist_of_rows);
}


void get_signatures_for_each_level(Sign signatures[], ConstantRowWeightSparseMatrix *matrix) {

    #if defined USE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS && USE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS == 0
        assert( matrix->nrows <= SBOX_LOG2_SIZE);
    #endif
    int arr[matrix->nrows];
    for (int i = 0; i < matrix->nrows; i++) {
        arr[i] = i;
        signatures[i] = get_signature_for_rows(matrix, arr, i + 1);
    }
}



// signatures[i][j] contains the signature of the matrix:
// [ L_W_A[0]   ]
// [ L_W_A[1]   ]
// [ ...        ]
// [ L_W_A[j-1] ]
// [ L_W_A[i]   ]
void get_all_signatures_for_each_level(Sign signatures[MAX_M][MAX_M], ConstantRowWeightSparseMatrix *matrix) {

    #if defined USE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS && USE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS == 0
        assert( matrix->nrows <= SBOX_LOG2_SIZE);
    #endif
    for (int i = 0; i < matrix->nrows; i++) {
        int arr[matrix->nrows];

        arr[0] = i;
        signatures[i][0] = get_signature_for_rows(matrix, arr, 1);

        for (int j = 1; j < matrix->nrows; j++) {
            arr[j - 1] = j - 1;
            arr[j] = i;
            signatures[i][j] = get_signature_for_rows(matrix, arr, j + 1);
            // printf("arr: \n");
            // for (int k = 0; k < j + 1; k++) {
                // printf("%d ", arr[k]);
            // }
            // printf("\n");
            // printf("signatures[%d][%d] = %u\n", i, j, signatures[i][j]);
        }
    }
}

void test_get_signatures_for_each_level(ConstantRowWeightSparseMatrix* matrix) {
    Sign signatures[matrix->nrows];
    get_signatures_for_each_level(signatures, matrix);
    for (int i = 0; i < matrix->nrows; i++) {
        printf("signatures[%d] = %x\n", i, signatures[i]);
    }
}


int is_a_valid_child(Sign signatures[], ConstantRowWeightSparseMatrix* target_matrix, int next_level, uint16_t base_sigs[], int candidate);

void find_all_possible_children(Sign signatures[], ConstantRowWeightSparseMatrix* target_matrix,
                                int current_support[], int next_level, int possibles[], int *ptr_npossible) {
    int iposs = 0;

    uint16_t sigs[target_matrix->ncols];
    memset(sigs, 0, target_matrix->ncols*sizeof(*sigs));

    get_signature_cummulative_base(target_matrix, current_support, next_level, sigs);

    for (int i = 0; i < *ptr_npossible; i++) {

        if (is_a_valid_child(signatures, target_matrix, next_level, sigs, i)) {
            possibles[iposs++] = possibles[i];
        }
    }
    *ptr_npossible = iposs;
    current_support[next_level] = -1;
}


int is_a_valid_child(Sign signatures[], ConstantRowWeightSparseMatrix* target_matrix, int next_level, uint16_t base_sigs[], int candidate) {
    uint16_t tmp_sigs[target_matrix->ncols];
    memcpy(tmp_sigs, base_sigs, target_matrix->ncols*sizeof(*base_sigs));

    Sign sig = get_signature_cummulative(target_matrix, candidate, next_level, tmp_sigs);
    return (signatures[next_level] == sig);
}


void attack_dfs_rec(PermutationExtractor* e, int level, PossibleChildrenPerLevel* poss_per_level);


int COUNT = 1;


int64_t SumTimes[100] = {0};
int64_t CountTimes[100] = {0};


int64_t PossiblesPerLevel[100] = {0};


void attack_dfs_rec(PermutationExtractor* e, int level, PossibleChildrenPerLevel* poss_per_level) {

    struct timespec start, end;
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &start);

    if (COUNT++ % 1000000 == 0) {
        double total = 1;
        int i;
        for (i = 0; CountTimes[i] == 0 && i < MAX_M; i++) {
            total *= PossiblesPerLevel[i];
        }
        printf("[*] Information on current level\n");
        printf("    Worst case total time: %lf seconds\n", total*SumTimes[i]/CountTimes[i]/1e9);
        printf("\n");
        printf("    Current path ");
        for (int i = 0; i < level; i++) {
            printf("%3d ", poss_per_level->support[i]);
        }
        printf("\n");
    }
    if (poss_per_level->level == poss_per_level->nlevels) {
        // This must run inside a critical section because M4RI cannot perform the
        // calculations with multiple threads. This can be solved by separating the memory
        // allocated for M4RI for each thread, but we did not do this in this version.
        #pragma omp critical
        {
            for (int i = 0; i < level; i++) {
                printf("%3d ", poss_per_level->support[i]);
            }
            printf("\n");
            printf("LEAF!\n");
            extract_permutation_from_matching(e, poss_per_level->target_matrix, poss_per_level->support, level);
            printf("    ... :^(\n");
        }
        return;
    }

    for (int ix = 0; ix < poss_per_level->nposs[level]; ix++) {
        int p = poss_per_level->poss[level][ix];

        poss_per_level_add_row_to_support(poss_per_level, p);
        int run = 1;
        for (int l = poss_per_level->level; l < poss_per_level->nlevels; l++) {
            if (poss_per_level->nposs[level] == 0) {
                run = 0;
                break;
            }
        }

        if (run) {
            attack_dfs_rec(e, level + 1, poss_per_level);
        }

        for (int l = level; l < poss_per_level->nlevels; l++) {
            poss_per_level->nposs[l] = poss_per_level->limits_of_removed_per_level[l][level];
        }
        poss_per_level->level--;
    }
    clock_gettime(CLOCK_THREAD_CPUTIME_ID, &end);

    int64_t delta_ns = ((end.tv_sec - start.tv_sec) * 1e9 + (end.tv_nsec - start.tv_nsec));

    SumTimes[level] += delta_ns;
    CountTimes[level]++;
}


// CAUTION: This function assumes the sparse entries are ordered from lowest to highest
int get_sparse_rows_hamming_distance(ConstantRowWeightSparseMatrix* target_matrix, int row_index1, int row_index2) {

    int equal_count = 0;
    for (int i = 0, j = 0; i < target_matrix->row_weight && j < target_matrix->row_weight;) {
        int v1 = target_matrix->rows[row_index1][i];
        int v2 = target_matrix->rows[row_index2][j];
        if (v1 == v2) {
            equal_count += 1;
            i += 1;
            j += 1;
        }

        i += (v1 < v2);
        j += (v2 < v1);
    }

    return 2*(target_matrix->row_weight - equal_count);
}


RowsPerDistance *precompute_rows_per_distance(ConstantRowWeightSparseMatrix* matrix, int row_index) {
    RowsPerDistance *rows_per_distance = malloc(sizeof(*rows_per_distance));

    rows_per_distance->n_for_distance = calloc(2*matrix->row_weight + 1, sizeof(*rows_per_distance->n_for_distance));

    for (int i = 0; i < matrix->nrows; i++) {
        int d = get_sparse_rows_hamming_distance(matrix, row_index, i);
        rows_per_distance->n_for_distance[d]++;
    }

    rows_per_distance->n = 2*matrix->row_weight + 1;
    rows_per_distance->data = malloc(rows_per_distance->n*sizeof(*rows_per_distance->data));
    for (int d = 0; d < rows_per_distance->n; d++) {
        rows_per_distance->data[d] = malloc((rows_per_distance->n_for_distance[d])*sizeof(*rows_per_distance->data[d]));
        rows_per_distance->n_for_distance[d] = 0;
    }

    for (int i = 0; i < matrix->nrows; i++) {
        int d = get_sparse_rows_hamming_distance(matrix, row_index, i);
        int j = rows_per_distance->n_for_distance[d]++;
        rows_per_distance->data[d][j] = i;
    }

    return rows_per_distance;
}

void free_rows_per_distance(RowsPerDistance *rows_per_distance) {
    for (int d = 0; d < rows_per_distance->n; d++)
        free(rows_per_distance->data[d]);
    free(rows_per_distance->n_for_distance);
    free(rows_per_distance);
}


PossibleChildrenPerLevel* init_poss_per_level(ConstantRowWeightSparseMatrix* target_matrix,
                                              ConstantRowWeightSparseMatrix* pattern_matrix) {
    PossibleChildrenPerLevel* poss_per_level;

    poss_per_level = malloc(sizeof(*poss_per_level));
    get_all_signatures_for_each_level(poss_per_level->signatures, pattern_matrix);

    // Sanity Check
    Sign x[pattern_matrix->nrows];
    get_signatures_for_each_level(x, pattern_matrix);

    for (int i = 0; i < pattern_matrix->nrows - 1; i++) {
        // printf("%d : %x -> %x\n", i, x[i], poss_per_level->signatures[i][i]);
        assert(x[i] == poss_per_level->signatures[i][i]);
    }

    poss_per_level->nlevels = pattern_matrix->nrows;
    for (int i = 0; i < poss_per_level->nlevels; i++) {
        for (int j = 0; j < poss_per_level->nlevels; j++) {
            poss_per_level->limits_of_removed_per_level[i][j] = target_matrix->nrows; // = nposs[l]
        }
    }

    poss_per_level->pattern_matrix = pattern_matrix;
    poss_per_level->target_matrix = target_matrix;
    poss_per_level->support = malloc(poss_per_level->nlevels*sizeof(*(poss_per_level->support)));
    poss_per_level->level = 0;
    poss_per_level->nposs = malloc(poss_per_level->nlevels*sizeof(*(poss_per_level->nposs)));
    for (int l = 0; l < poss_per_level->nlevels; l++) {
        poss_per_level->nposs[l] = 0;
    }
    poss_per_level->poss = malloc(poss_per_level->nlevels*sizeof(*(poss_per_level->poss)));
    for (int l = 0; l < poss_per_level->nlevels; l++) {
        poss_per_level->nposs[l] = target_matrix->nrows;
        poss_per_level->poss[l] = malloc(poss_per_level->nposs[l]*sizeof(*(poss_per_level->poss[l])));
        for (int j = 0; j < poss_per_level->nposs[l]; j++) {
            poss_per_level->poss[l][j] = j;
        }
    }

    return poss_per_level;
}


void poss_per_level_add_row_to_support(PossibleChildrenPerLevel* poss_per_level, int row) {
    assert(row < poss_per_level->target_matrix->nrows);
    assert(poss_per_level->level < poss_per_level->nlevels);

    poss_per_level->support[poss_per_level->level] = row;
    poss_per_level->level++;

    uint16_t base_sigs[poss_per_level->target_matrix->ncols];
    memset(base_sigs, 0, poss_per_level->target_matrix->ncols*sizeof(*base_sigs));
    get_signature_cummulative_base(poss_per_level->target_matrix, poss_per_level->support, poss_per_level->level, base_sigs);

    for (int level = poss_per_level->level; level < poss_per_level->nlevels; level++) {
        int base_pos = poss_per_level->limits_of_removed_per_level[level][poss_per_level->level - 1];
        poss_per_level->limits_of_removed_per_level[level][poss_per_level->level] = base_pos;

        for (int i = 0; i < poss_per_level->nposs[level];) {
            if (poss_per_level->limits_of_removed_per_level[level][poss_per_level->level] == 0) {
                poss_per_level->nposs[level] = 0;
                break;
            }

            uint16_t tmp_sigs[poss_per_level->target_matrix->ncols];
            memcpy(tmp_sigs, base_sigs, poss_per_level->target_matrix->ncols*sizeof(*base_sigs));

            int p = poss_per_level->poss[level][i];

            Sign sig = get_signature_cummulative(poss_per_level->target_matrix, p, poss_per_level->level, tmp_sigs);

            // Should remove p putting it in the end
            if (poss_per_level->signatures[level][poss_per_level->level] != sig) {

                poss_per_level->limits_of_removed_per_level[level][poss_per_level->level]--;

                int target = poss_per_level->limits_of_removed_per_level[level][poss_per_level->level];
                int tmp = poss_per_level->poss[level][target];
                poss_per_level->poss[level][target] = poss_per_level->poss[level][i];
                poss_per_level->poss[level][i] = tmp;

                poss_per_level->nposs[level]--;
            }
            else {
                i++;
            }
        }
    }

    // printf("=======================================\n");
    // for (int l = 0; l < poss_per_level->nlevels; l++) {
    //     printf("LEVEL %d: [%d] -> %d\n", poss_per_level->level, l,
    //             poss_per_level->limits_of_removed_per_level[l][poss_per_level->level]);
    // }
}

void attack_dfs(PermutationExtractor* e,
                ConstantRowWeightSparseMatrix* pattern_matrix,
                ConstantRowWeightSparseMatrix* target_matrix,
                int secret_support[]) {


    Sign signatures[MAX_M][MAX_M];
    get_all_signatures_for_each_level(signatures, pattern_matrix);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < target_matrix->nrows; i++) {
        PossibleChildrenPerLevel* poss_per_level = init_poss_per_level(target_matrix, pattern_matrix);
        poss_per_level_add_row_to_support(poss_per_level, i);
        attack_dfs_rec(e, 1, poss_per_level);
    }
}

void remove_from_vector(int vector[], int* ptr_nvector, int index) {
    vector[index] = vector[--(*ptr_nvector)];
}


void get_npossibles_per_level(ConstantRowWeightSparseMatrix* pattern_matrix,
                              ConstantRowWeightSparseMatrix* target_matrix,
                              int secret_support[]) {
    int current_support[pattern_matrix->nrows];
    int npossibles[pattern_matrix->nrows];

    Sign signatures[pattern_matrix->nrows];

    get_signatures_for_each_level(signatures, pattern_matrix);

    /* sanity check */
    for (int i = 0; i < pattern_matrix->nrows; i++) {
        current_support[i] = secret_support[i];
        Sign sig = get_signature_for_rows(target_matrix, current_support, i + 1);
        assert(signatures[i] == sig);
    }

    npossibles[0] = target_matrix->nrows;
    for (int level = 1; level < pattern_matrix->nrows; level++) {
        int npossible = target_matrix->nrows;
        int possibles[target_matrix->nrows];

        for (int i = 0; i < target_matrix->nrows; i++)
            possibles[i] = i;

        current_support[level - 1] = secret_support[level - 1];
        find_all_possible_children(signatures, target_matrix, current_support, level, possibles, &npossible);
        npossibles[level] = npossible;
        // printf("npossibles[%d] = %d\n", level, npossibles[level]);

    }

    for (int i = 0; i < pattern_matrix->nrows; i++) {
        printf("npossibles[%d] = %d\n", i, npossibles[i]);
        PossiblesPerLevel[i] = npossibles[i];
    }
}


int test_main(int argc, char* argv[]) {

    int ma = 10;
    ConstantRowWeightSparseMatrix* matrix = read_sparse_matrix_input(argv[1], ma);
    test_get_signatures_for_each_level(matrix);
    test_read_sparse_matrix_input(argv[1]);

    return 0;
}


Label *compute_col_labels(mzd_t* matrix) {
    Label *labels = calloc(matrix->ncols, sizeof(*labels));
    for (int j = 0; j < matrix->ncols; j++) {
        for (int i = 0; i < matrix->nrows; i++) {
            if (mzd_read_bit(matrix, i, j))
                labels[j] += 1 << i;
        }
        printf("labels[%d] = %d\n", j, labels[j]);
    }
    return labels;
}

void compute_labels_to_indexes(PermutationExtractor *pe) {
    for (int i = 0; i < MAX_LABEL; i++) {
        pe->labels_to_Aindexes[i].n = 0;
    }

    for (int i = 0; i < pe->A->ncols; i++) {
        pe->labels_to_Aindexes[pe->Alabels[i]].n++;
    }

    for (int i = 0; i < MAX_LABEL; i++) {
        int n = pe->labels_to_Aindexes[i].n;
        pe->labels_to_Aindexes[i].values = malloc(n*sizeof(*(pe->labels_to_Aindexes[i].values)));
        pe->labels_to_Aindexes[i].n = 0;
    }

    for (int i = 0; i < pe->A->ncols; i++) {
        int n = pe->labels_to_Aindexes[pe->Alabels[i]].n;
        pe->labels_to_Aindexes[pe->Alabels[i]].values[n] = i;
        pe->labels_to_Aindexes[pe->Alabels[i]].n++;
    }

    // #if DEBUG
    for (int i = 0; i < MAX_LABEL; i++) {
        if (pe->labels_to_Aindexes[i].n > 0) {
            printf("%d: ", i);
            for (int j = 0; j < pe->labels_to_Aindexes[i].n; j++) {
                printf("%d ", pe->labels_to_Aindexes[i].values[j]);
            }
            printf("\n");
        }
    }
    // #endif
}

typedef struct _int_pair {
    int index;
    int r;
} IntPair;


static int int_pair_cmp(const void *p1, const void *p2) {
    return (*(IntPair *)p1).r - (*(IntPair *)p1).r;
}


void compute_A1_A2_indexes_and_repetitions(PermutationExtractor* pe) {
    IntPair ind_rs[pe->A->ncols];
    int k = 0;
    for (int i = 0; i < MAX_LABEL; i++) {
        for (int j = 0; j < pe->labels_to_Aindexes[i].n; j++) {
            ind_rs[k].index = pe->labels_to_Aindexes[i].values[j];
            ind_rs[k].r = pe->labels_to_Aindexes[i].n;
            k++;
        }
    }
    qsort(ind_rs, pe->A->ncols, sizeof(*ind_rs), int_pair_cmp);
    #if DEBUG
        for (int i = 0; i < pe->A->ncols; i++) {
            printf("ind:%d r:%d\n", ind_rs[i].index, ind_rs[i].r);
        }
        mzd_print(pe->A);
    #endif
    printf("=======================\n");
    mzd_t *sorted_matrix = mzd_init(pe->A->nrows, pe->A->ncols);
    for (int i = 0; i < pe->A->nrows; i++) {
        for (int j = 0; j < pe->A->ncols; j++) {
            int perm_bit = mzd_read_bit(pe->A, i, ind_rs[j].index);
            mzd_write_bit(sorted_matrix, i, j, perm_bit);
        }
    }

    mzd_echelonize(sorted_matrix, 0);


    pe->A2_indexes = malloc(pe->A->nrows*sizeof(*pe->A2_indexes));
    pe->A1_indexes = malloc((pe->A->ncols - pe->A->nrows)*sizeof(*pe->A1_indexes));

    int na1 = 0;
    int na2 = 0;
    int i = 0;
    for (int j = 0; j < pe->A->ncols; j++) {
        if (i < pe->A->nrows && mzd_read_bit(sorted_matrix, i, j) == 1) {
            pe->A2_indexes[na2++] = ind_rs[j].index;
            i++;
            printf("pe->A2_indexes[%d] = %d\n", na2 - 1, pe->A2_indexes[na2 - 1]);
        }
        else {
            pe->A1_indexes[na1++] = ind_rs[j].index;
            printf("pe->A1_indexes[%d] = %d\n", na1 - 1, pe->A1_indexes[na1 - 1]);
        }
    }
    pe->n_a1 = na1;
    pe->n_a2 = na2;

    for (int l = 0; l < MAX_LABEL; l++)
        pe->repetitions[l] = 0;

    for (int i = 0; i < na1; i++) {
        int l = pe->Alabels[pe->A1_indexes[i]];
        pe->repetitions[l]++;
    }

    /* Start from labels_to_Aindexes */
    for (int l = 0; l < MAX_LABEL; l++) {
        pe->labels_to_A1indexes[l].n = 0;
        int n = pe->labels_to_Aindexes[l].n;
        pe->labels_to_A1indexes[l].values = malloc(n*sizeof(*pe->labels_to_A1indexes[l].values));
    }
    for (int i = 0; i < na1; i++) {
        int l = pe->Alabels[pe->A1_indexes[i]];
        int n = pe->labels_to_A1indexes[l].n;
        pe->labels_to_A1indexes[l].values[n] = i;
        pe->labels_to_A1indexes[l].n = n + 1;
    }

    #if DEBUG
        for (int i = 0; i < na1; i++) {
            int l = pe->Alabels[pe->A1_indexes[i]];
            printf("A1[%d] = %d -> group r = %d\n", i, pe->A1_indexes[i], pe->repetitions[l]);
        }

        mzd_print(pe->A);
        printf("====================\n");
        mzd_print(sorted_matrix);
    #endif

    mzd_free(sorted_matrix);
}

void build_A2_inv_times_A1(PermutationExtractor* pe) {
    mzd_t* A2 = mzd_init(pe->n_a2, pe->n_a2);
    for (int i = 0; i < pe->n_a2; i++) {
        for (int j = 0; j < pe->n_a2; j++) {
            int bit = mzd_read_bit(pe->A, i, pe->A2_indexes[j]);
            mzd_write_bit(A2, i, j, bit);
        }
    }
    mzd_t* A2_inv = mzd_init(pe->n_a2, pe->n_a2);
    mzd_inv_m4ri(A2_inv, A2, 0);

    printf("Matrix A1\n");
    mzd_print(pe->A1);

    pe->A2_inv_times_A1 = mzd_mul(NULL, A2_inv, pe->A1, 0);

    #if DEBUG
        printf("A2inv\n");
        mzd_print(pe->A2_inv_times_A1);
    #endif

    mzd_free(A2);
    mzd_free(A2_inv);

}

void compute_unique_labels(PermutationExtractor* pe) {
    int n = 0;
    for (int i = 0; i < MAX_LABEL; i++) {
        if (pe->labels_to_A1indexes[i].n > 0)
            n++;
    }
    pe->unique_labels = malloc(n*sizeof(*pe->unique_labels));
    n = 0;
    for (int i = 0; i < MAX_LABEL; i++) {
        if (pe->labels_to_A1indexes[i].n > 0)
            pe->unique_labels[n++] = i;
    }
    pe->n_unique_labels = n;
}

void build_A1(PermutationExtractor* pe) {
    pe->A1 = mzd_init(pe->A->nrows, pe->n_a1);
    for (int i = 0; i < pe->A->nrows; i++) {
        for (int j = 0; j < pe->n_a1; j++) {
            int bit = mzd_read_bit(pe->A, i, pe->A1_indexes[j]);
            mzd_write_bit(pe->A1, i, j, bit);
        }
    }
    #if DEBUG
        printf("A1\n");
        mzd_print(pe->A1);
    #endif
}

void compute_Vpub_signatures(PermutationExtractor* pe) {
    for (int i = 0; i < pe->Vpub->nrows; i++) {
        pe->Vpub_signatures[i] = 0;
        for (int j = 0; j < pe->Vpub->ncols; j++) {
            pe->Vpub_signatures[i] += (mzd_read_bit(pe->Vpub, i, j) << j);
        }
    }
    for (int i = 0; i < MAX_VPUB_SIGNATURES; i++) {
        pe->Vpub_signatures_count[i] = 0;
    }
    for (int i = 0; i < pe->Vpub->nrows; i++) {
        int signature = pe->Vpub_signatures[i];
        pe->Vpub_signatures_count[signature]++;
    }
}

// CAUTION: We assume that matrix A is in the form
//
//      A = [L_A]
//          [A' ]
//  Where L_A is the set of n_A_low L.I. low weight keywords of A.
//  This may be changed in the future.
PermutationExtractor* init_permutation_extractor(char A_jcf_file[], char Vpub_jcf_file[], int n_A_low) {
    PermutationExtractor* pe = malloc(sizeof(*pe));

    pe->A = mzd_from_jcf(A_jcf_file, 1);
    pe->Vpub = mzd_from_jcf(Vpub_jcf_file, 1);

    pe->Alow = mzd_submatrix(NULL, pe->A, 0, 0, n_A_low, pe->A->ncols);
    printf("A_low_weight_head\n");
    mzd_print(pe->A);
    printf("A_LWSA\n");
    mzd_print(pe->Alow);
    printf("Vpub\n");
    mzd_print(pe->Vpub);

    pe->Alabels = compute_col_labels(pe->Alow);
    compute_labels_to_indexes(pe);
    compute_A1_A2_indexes_and_repetitions(pe);
    build_A1(pe);
    build_A2_inv_times_A1(pe);
    compute_unique_labels(pe);
    compute_Vpub_signatures(pe);

    printf("======================= %d ---\n", pe->n_a2);
    for (int i = 0; i < pe->n_a2; i++) {
        printf("pe->A2_indexes[%d] = %d (%d)\n", i, pe->A2_indexes[i], pe->labels_to_Aindexes[pe->Alabels[pe->A2_indexes[i]]].n);
    }
    printf("NEXT\n");
    fflush(stdout);
    for (int i = 0; i < pe->n_a1; i++) {
        printf("pe->A1_indexes[%d] = %d (%d)\n", i, pe->A1_indexes[i], pe->labels_to_Aindexes[pe->Alabels[pe->A1_indexes[i]]].n);
        fflush(stdout);

    }

    return pe;
}


Label* compute_col_labels_sparse(ConstantRowWeightSparseMatrix* M, int* rowset, int n_rowset) {
    Label *labels = calloc(M->ncols, sizeof(*labels));
    for (int i = 0; i < n_rowset; i++) {
        int row = rowset[i];
        for (int j = 0; j < M->row_weight; j++) {
            labels[M->rows[row][j]] += (1 << i);
        }
    }
    return labels;
}

void extract_permutation_from_matching_rec(PermutationExtractor* e, List labels_to_Vpub_indexes[], mzd_t* Vsec_partial, int label_index, int start, int nrows_set);

void extract_permutation_from_matching(PermutationExtractor* e, ConstantRowWeightSparseMatrix* M, int* rowset, int n_rowset) {
    Label* target_labels;
    // for (int i = 0; i < e->n_a1; i++) {
    //     Label l = e->Alabels[e->A1_indexes[i]];
    //     printf("%d -> %d\n", l, e->labels_to_Aindexes[l].n);
    // }

    target_labels = compute_col_labels_sparse(M, rowset, n_rowset);
    // for (int i = 0; i < M->ncols; i++) {
    //     printf("target labels[%d] = %d\n", i, target_labels[i]);
    // }

    mzd_t* Vsec_partial = mzd_init(e->n_a1, e->Vpub->ncols);
    int nused[MAX_LABEL] = {0};
    List labels_to_Vpub_indexes[MAX_LABEL];

    for (int i = 0; i < e->Vpub->nrows; i++) {
        Label row_label = target_labels[i];
        nused[row_label]++;
    }

    for (int i = 0; i < MAX_LABEL; i++) {
        int n = nused[i];
        labels_to_Vpub_indexes[i].values = malloc(n*sizeof(*labels_to_Vpub_indexes[i].values));
        labels_to_Vpub_indexes[i].n = n;
        nused[i] = 0;
    }
    for (int i = 0; i < e->Vpub->nrows; i++) {
        Label row_label = target_labels[i];
        labels_to_Vpub_indexes[row_label].values[nused[row_label]] = i;
        nused[row_label]++;
    }
    // mzd_print(Vsec_partial);
    // printf("==================\n");
    // mzd_print(e->Vpub);
    // exit(1);
    extract_permutation_from_matching_rec(e, labels_to_Vpub_indexes, Vsec_partial, 0, 0, 0);
    mzd_free(Vsec_partial);
}


void swap(int *a, int *b) {
    int tmp;

    tmp = *a;
    *a = *b;
    *b = tmp;
}


int is_possible_key(PermutationExtractor* e, mzd_t* Vsec_partial, mzd_t* product) {
    int candidate_signatures_count[MAX_VPUB_SIGNATURES] = {0};

    for (int i = 0; i < product->nrows; i++) {
        int signature = 0;
        for (int j = 0; j < product->ncols; j++)
            signature += mzd_read_bit(product, i, j) << j;
        if (e->Vpub_signatures_count[signature] == 0)
            return 0;

        candidate_signatures_count[signature]++;
    }

    for (int i = 0; i < Vsec_partial->nrows; i++) {
        int signature = 0;
        for (int j = 0; j < Vsec_partial->ncols; j++)
            signature += mzd_read_bit(Vsec_partial, i, j) << j;

        candidate_signatures_count[signature]++;
    }
    for (int i = 0; i < MAX_VPUB_SIGNATURES; i++) {
        if (candidate_signatures_count[i] != e->Vpub_signatures_count[i])
            return 0;
    }
    return 1;
}


void extract_permutation_from_matching_rec(PermutationExtractor* e, List labels_to_Vpub_indexes[], mzd_t* Vsec_partial, int label_index, int start, int nrows_set) {
    if (nrows_set >= e->n_a1) {
        for (int li = 0; li < label_index; ++li) {
            int label = e->unique_labels[li];
            for (int i = 0; i < e->labels_to_A1indexes[label].n; i++) {
                mzd_copy_row(Vsec_partial, e->labels_to_A1indexes[label].values[i],
                             e->Vpub, labels_to_Vpub_indexes[label].values[i]);
            }
        }

        mzd_t* product = mzd_mul(NULL, e->A2_inv_times_A1, Vsec_partial, 0);

        if (is_possible_key(e, Vsec_partial, product)) {
            printf("---- PARTIALS --------\n");
            mzd_print(Vsec_partial);
            printf("-------------\n");
            mzd_print(product);
            printf("---- END PARTIALS --------\n");



            printf("================== Vsec = KEY ===========\n");

            mzd_t* Vsec = mzd_init(e->Vpub->nrows, e->Vpub->ncols);

            for (int i = 0; i < e->n_a2; i++) {
                mzd_copy_row(Vsec, e->A2_indexes[i],
                             product, i);
            }
            for (int i = 0; i < e->n_a1; i++) {
                mzd_copy_row(Vsec, e->A1_indexes[i],
                             Vsec_partial, i);
            }
            mzd_print(Vsec);
            mzd_free(Vsec);

            printf("============== KEY FOUND ====================\n");
            exit(1);
        }
        mzd_free(product);
        return;
    }
    int label = e->unique_labels[label_index];

    if (start == e->repetitions[label]) {
        extract_permutation_from_matching_rec(e, labels_to_Vpub_indexes, Vsec_partial, label_index + 1, 0, nrows_set + e->repetitions[label]);
    }

    int *arr = labels_to_Vpub_indexes[label].values;
    for (int i = start; i < labels_to_Vpub_indexes[label].n; i++) {
        swap(arr + i, arr + start);
        extract_permutation_from_matching_rec(e, labels_to_Vpub_indexes, Vsec_partial, label_index, start + 1, nrows_set);
        swap(arr + i, arr + start);
    }
}


// If you are using Windows, make sure to change the directory separator below
void get_full_path_for_file(char output[], char directory[], char filename[]) {
    strcpy(output, directory);
    strcat(output, "/" );
    strcat(output, filename);
}


int read_param_la(char filepath[]) {
    int la;
    FILE *file = fopen(filepath, "r");
    if (!file) {
        fprintf(stderr, "Could not open param_la file: %s.\n", filepath);
        exit(1);
    }
    fscanf(file, "%d", &la);
    fclose(file);

    return la;
}


int main(int argc, char* argv[]) {

    if (argc != 2) {
        printf("Usage: %s <path to directory containing challenge files>\n", argv[0]);
        return 1;
    }

    char filepaths[6][1000];

    get_full_path_for_file(filepaths[0], argv[1], "la.param");
    get_full_path_for_file(filepaths[1], argv[1], "LWSA.sparse_matrix");
    get_full_path_for_file(filepaths[2], argv[1], "LWSK.sparse_matrix");
    get_full_path_for_file(filepaths[3], argv[1], "secret_subset.indexes");
    get_full_path_for_file(filepaths[4], argv[1], "A_low_weight_head.jcf");
    get_full_path_for_file(filepaths[5], argv[1], "Vpub.jcf");

    int ma = read_param_la(filepaths[0]);


    ConstantRowWeightSparseMatrix* pattern_matrix = read_sparse_matrix_input(filepaths[1], ma);
    ConstantRowWeightSparseMatrix* target_matrix = read_sparse_matrix_input(filepaths[2], -1);
    int secret_support[ma];
    read_secret_support_file(secret_support, ma, filepaths[3]);

    PermutationExtractor* e = init_permutation_extractor(filepaths[4], filepaths[5], ma);
    printf("Each permutation extraction will take:\n");
    for (int i = 0; i < e->n_unique_labels; i++) {
        int n = e->labels_to_Aindexes[e->unique_labels[i]].n;
        if (n > 1) {
            printf("(%d!/%d!)", n, n - e->repetitions[e->unique_labels[i]]);
        }
    }
    printf("\n");

    get_npossibles_per_level(pattern_matrix, target_matrix, secret_support);

    // To perform a sanity check, you can uncomment the following line of code, which will
    // try to extract the secret permutation of the correct matching. Notice that if you
    // are using attack_param_la and attack_param_w too small, this will take a LIFETIME
    // TO COMPLETE!
    // extract_permutation_from_matching(e, target_matrix, secret_support, ma);

    attack_dfs(e, pattern_matrix, target_matrix, secret_support);
    free(pattern_matrix);
    free(target_matrix);
    return 0;
}
