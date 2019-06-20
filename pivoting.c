/********************************************************************
 *
 * pivoting.c
 *
 * This library provides routines to pivot matrices of a linear
 * system of equations. Solving the system would require much lower
 * computational cost after pivoting. Two pivoting algorithms, namely
 * inactivation and Zlatev pivoting, are implemented. One interface
 * is exposed to caller.
 *
 ********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include "galois.h"
#define ZLATEVS  3    // Number of searching rows using Zlatev's strategy
#define IA_INIT  0    // Initial number of inactivated columns
#define IA_STEP  1    // Gradually inactivate more columns

#define TRACE    5    // trace log level
static int loglevel = 5;    // log level for the library
int get_loglevel()
{
    return loglevel;
}
/*
 * Pivoting algorithms use double-linked lists to store numbers of
 * nonzeros of rows and columns of a matrix.
 */
typedef struct subscript  Subscript;
typedef struct subscripts ssList;
struct subscript
{
    int index;
    int nonzeros;
    struct subscript *next;
    struct subscript *prev;
};

struct subscripts
{
    struct subscript *ssFirst;
    struct subscript *ssLast;
};



/*
 * Procedures to pivot matrix.
 */
static int inactivation_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots);
static int zlatev_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots);

/*
 * Reshape matrix A and B of Ax=B according to pivot sequence in (RowPivots, ColPivots)
 */
static void reshape_matrix(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, ssList *RowPivots, ssList *ColPivots);

/*
 * Reorder columns of matrix according to sequence in ColPivots
 */
static void permute_matrix_columns(int nrow, int ncolA, GF_ELEMENT **A, ssList *ColPivots);

/*
 * Helper functions for pivoting
 * We use double-linked lists when pivoting a matrix.
 */
static void insertSubAtBeginning(ssList *sub_list, Subscript *sub);
static void insertSubAtEnd(ssList *sub_list, Subscript *sub);
static void removeSubscript(ssList *sub_list, Subscript *sub);
static void free_subscriptList(ssList *sub_list);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);

/**********************************************************************************
 * pivot_matrix_x()
 * Main interfaces to the library.
 * Input:
 *  Ax = B
 *  It is required that 2-D matrices A and B are stored in form of double pointers
 *
 *  A ->A[0] A[0][1] A[0][2] ...
 *      A[1] ...
 *      .
 *      .
 * Parameters:
 *  nrow  - number of rows of A and B
 *  ncolA - number of columns of A
 *  ncolB - number of columns of B
 *
 * Return:
 *  number of Galois field operations consumed
 *
 * Return as arguments:
 *  ctoo_r/ctoo_c  - Arrays containing mappings of row/col indices after pivoting
 *  inactives      - number of inactivated columns
 *
 *  [A] will be transformed into the form of:
 *      -                 -
 *      | x 0 0 0 0 x x x |
 *      | 0 x 0 0 0 x x x |
 *      | 0 0 x 0 0 x x x |
 *      | 0 0 0 x 0 x x x |
 *      | 0 0 0 0 x x x x |
 *      | - - - - - ------|
 *      | 0 0 0 0 0|x x x |
 *      | 0 0 0 0 0|0 x x |
 *      | 0 0 0 0 0|0 0 x |
 *      -                 -
 *  and B will be processed accordingly.
 *
 **********************************************************************************/

long pivot_matrix_oneround(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int **ptr_ctoo_r, int **ptr_ctoo_c, int *inactives)
{
    int i, j, k;
    long operations = 0;
    ssList *row_pivots = malloc(sizeof(ssList));
    ssList *col_pivots = malloc(sizeof(ssList));
    row_pivots->ssFirst = row_pivots->ssLast = NULL;
    col_pivots->ssFirst = col_pivots->ssLast = NULL;

    double pivoting_time=0.0;
    clock_t start_pivoting, stop_pivoting;
    long long nonzeros, total;
    if (get_loglevel() == TRACE) {
        printf("Start one round of pivoting...\n");
        printf("Inactivation pivoting... IA_INIT: %d, IA_STEP: %d.\n", IA_INIT, IA_STEP);
        start_pivoting = clock();
    }
    int ias = inactivation_pivoting(nrow, ncolA, A, row_pivots, col_pivots);
    *inactives = ias;
    if (get_loglevel() == TRACE) {
        printf("A total of %d/%d columns are inactivated.\n", ias, ncolA);
        stop_pivoting = clock();
        pivoting_time = ((double) (stop_pivoting - start_pivoting)) / CLOCKS_PER_SEC;
        printf("Inactivation pivoting consumed time: %.6f seconds.\n", pivoting_time);
    }

    // Save orders of row/column indices after pivoting
    // Original-to-Current mapping: the i-th element of the array gives where is the original i-th row/col in the new (virtually) re-ordered matrix
    // Current-to-Original mapping: the i-th element of the array specifies what was the original row/col id of the current i-th row/col
    int *ctoo_r = *ptr_ctoo_r;
    int *ctoo_c = *ptr_ctoo_c;
    Subscript *ss_pt_r = row_pivots->ssFirst;
    Subscript *ss_pt_c = col_pivots->ssFirst;
    for (i=0; i<ncolA; i++) {
        ctoo_r[i] = ss_pt_r->index;
        ss_pt_r = ss_pt_r->next;
        ctoo_c[i] = ss_pt_c->index;
        ss_pt_c = ss_pt_c->next;
    }
    free_subscriptList(row_pivots);
    free_subscriptList(col_pivots);

    // Diagonalize active part
    // Don't physically swap row/col. Perform all operations on the original matrix with the help of the mappings 
    // To save random access time, we still need to make a copy of the top-right (nColA-ias) x ias matrix U.
    GF_ELEMENT **U = calloc(ncolA-ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ncolA-ias; i++) {
        U[i] = calloc(ias, sizeof(GF_ELEMENT));
        for (j=0; j<ias; j++)
            U[i][j] = A[ctoo_r[i]][ctoo_c[ncolA-ias+j]];
    }
    // And also make a copy of the bottom-right (ias x ias) matrix T.
    // This kind of copy is a tradeoff between memory usage and access time. 
    GF_ELEMENT **T = calloc(ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ias; i++){
        T[i] = calloc(ias, sizeof(GF_ELEMENT));
        for (j=0; j<ias; j++)
            T[i][j] = A[ctoo_r[ncolA-ias+i]][ctoo_c[ncolA-ias+j]];
    }
    // Start to diagonalize
    clock_t start_p, stop_p;
    start_p = clock();
    long long ops1=0;
    GF_ELEMENT quotient;
    nonzeros = 0;
    for (i=0; i<ncolA-ias; i++) {
        for (j=i+1; j<ncolA; j++) {
            if (A[ctoo_r[i]][ctoo_c[i]] == 0 && get_loglevel() == TRACE)
                printf("The diagonal element after re-ordering is nonzero.\n");
            // process the item on (j, i)
            if (A[ctoo_r[j]][ctoo_c[i]] != 0) {
                quotient = galois_divide(A[ctoo_r[j]][ctoo_c[i]], A[ctoo_r[i]][ctoo_c[i]]);
                ops1 += 1;
                // multiply-and-add the corresponding part in the inactive part
                if (j < ncolA-ias) {
                    // eliminating nonzeros in the first ncolA-ias rows (but below diagonal)
                    galois_multiply_add_region(U[j], U[i], quotient, ias);
                } else {
                    // eliminating nonzeros in the last ias rows
                    galois_multiply_add_region(T[j-(ncolA-ias)], U[i], quotient, ias);
                }
                ops1 += ias;    // This part of matrix in processing is lower triangular in part, so operations only needed in the back half (i.e., inactiavted part)
                // simultaneously do the same thing on right matrix B
                galois_multiply_add_region(B[ctoo_r[j]], B[ctoo_r[i]], quotient, ncolB);
                ops1 += ncolB;
                A[ctoo_r[j]][ctoo_c[i]] = 0;            // eliminate the item
                nonzeros += 1;
            }
        }
    }
    stop_p = clock();
    if (get_loglevel() == TRACE) {
        printf("Diagonalize active part (%lld nonzero elements) after inactivating %d took %.6f seconds, cost %lld operations.\n", nonzeros, ias, ((double) (stop_p-start_p))/CLOCKS_PER_SEC, ops1);
    }
    operations += ops1;

    // Copy U and T back, free U but keep T for forward substitution
    for (i=0; i<ncolA-ias; i++) {
        for (j=0; j<ias; j++)
            A[ctoo_r[i]][ctoo_c[ncolA-ias+j]] = U[i][j];
        free(U[i]);
    }
    free(U);
    for (i=0; i<ias; i++) {
        for (j=0; j<ias; j++)
            A[ctoo_r[ncolA-ias+i]][ctoo_c[ncolA-ias+j]] = T[i][j];
    }

    /* Perform forward substitution on the ias x ias dense inactivated matrix. */
    start_p = clock();
    /*
    if (get_loglevel() == TRACE) {
        nonzeros = 0;
        total = 0;
        for (i=0; i<ias; i++) {
            for (j=0; j<ias; j++) {
                total++;
                if (T[i][j] != 0)
                    nonzeros++;
            }
        }
        printf("Current percentage of nonzeros of the sub-matrix: %.2f%%\n", (double) nonzeros/total*100);
        printf("Perform forward substitution on the square sub-matrix of size: %d x %d\n", ias, ias);
    }
    */
    // Make a copy of the corresponding msg matrices of T before performing forward substitution. 
    GF_ELEMENT **msg_submatrix = calloc(ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ias; i++){
        msg_submatrix[i] = calloc(ncolB, sizeof(GF_ELEMENT));
        memcpy(msg_submatrix[i], B[ctoo_r[ncolA-ias+i]], ncolB*sizeof(GF_ELEMENT));
    }

    long long ops = forward_substitute(ias, ias, ncolB, T, msg_submatrix);
    operations += ops;
    // Save the processed inactivated part back to A
    for (i=0; i<ias; i++) {
        for (j=0; j<ias; j++)
            A[ctoo_r[ncolA-ias+i]][ctoo_c[ncolA-ias+j]] = T[i][j];
        memcpy(B[ctoo_r[ncolA-ias+i]], msg_submatrix[i], ncolB*sizeof(GF_ELEMENT));
    }

    //free T, msg_submatrix
    for (i=0; i<ias; i++) {
        free(T[i]);
        free(msg_submatrix[i]);
    }
    free(T);
    free(msg_submatrix);
    stop_p = clock();
    if (get_loglevel() == TRACE) {
        printf("Forward substitution on the bottom-right part took %.6f seconds, cost %lld operations\n", ((double) (stop_p-start_p))/CLOCKS_PER_SEC, ops);
    }

    /*
    if (get_loglevel() == TRACE) {
        int missing_pivots = 0;
        for (i=0; i<ncolA; i++) {
            if (A[ctoo_r[i]][ctoo_c[i]] == 0)
                missing_pivots += 1;
        }
        printf("There are %d pivots missing after forward substitution\n", missing_pivots);
    }
    */

    return operations;
}

/********************************************************************************
 *
 * Pivot matrix for two rounds. First round uses inactivation pivoting against
 * the whole matrix. The second round performs on the relatively denser
 * bottom-right corner corresponding to the inactivated columns. Zlatev pivoting
 * is used for the second round.
 *
 ********************************************************************************/

long pivot_matrix_tworound(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int **ptr_ctoo_r, int **ptr_ctoo_c, int *inactives)
{
    int i, j, k;
    long operations = 0;
    // First pivoting: inactivation
    ssList *row_pivots = malloc(sizeof(ssList));
    ssList *col_pivots = malloc(sizeof(ssList));
    row_pivots->ssFirst = row_pivots->ssLast = NULL;
    col_pivots->ssFirst = col_pivots->ssLast = NULL;

    double pivoting_time=0.0;
    clock_t start_pivoting, stop_pivoting;
    long long nonzeros, total;
    if (get_loglevel() == TRACE) {
        printf("Start one round of pivoting...\n");
        printf("Inactivation pivoting... IA_INIT: %d, IA_STEP: %d.\n", IA_INIT, IA_STEP);
        start_pivoting = clock();
    }
    int ias = inactivation_pivoting(nrow, ncolA, A, row_pivots, col_pivots);
    *inactives = ias;
    if (get_loglevel() == TRACE) {
        printf("A total of %d/%d columns are inactivated.\n", ias, ncolA);
        stop_pivoting = clock();
        pivoting_time = ((double) (stop_pivoting - start_pivoting)) / CLOCKS_PER_SEC;
        printf("Inactivation pivoting consumed time: %.6f seconds.\n", pivoting_time);
    }

    // Save orders of row/column indices after pivoting
    // Original-to-Current mapping: the i-th element of the array gives where is the original i-th row/col in the new (virtually) re-ordered matrix
    // Current-to-Original mapping: the i-th element of the array specifies what was the original row/col id of the current i-th row/col
    int *ctoo_r = *ptr_ctoo_r;
    int *ctoo_c = *ptr_ctoo_c;
    Subscript *ss_pt_r = row_pivots->ssFirst;
    Subscript *ss_pt_c = col_pivots->ssFirst;
    for (i=0; i<ncolA; i++) {
        ctoo_r[i] = ss_pt_r->index;
        ss_pt_r = ss_pt_r->next;
        ctoo_c[i] = ss_pt_c->index;
        ss_pt_c = ss_pt_c->next;
    }
    free_subscriptList(row_pivots);
    free_subscriptList(col_pivots);

    // Diagonalize active part
    // Don't physically swap row/col. Perform all operations on the original matrix with the help of the mappings 
    // To save random access time, we still need to make a copy of the top-right (nColA-ias) x ias matrix U.
    GF_ELEMENT **U = calloc(ncolA-ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ncolA-ias; i++) {
        U[i] = calloc(ias, sizeof(GF_ELEMENT));
        for (j=0; j<ias; j++)
            U[i][j] = A[ctoo_r[i]][ctoo_c[ncolA-ias+j]];
    }
    // And also make a copy of the bottom-right (ias x ias) matrix T.
    // This kind of copy is a tradeoff between memory usage and access time. 
    GF_ELEMENT **T = calloc(ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ias; i++){
        T[i] = calloc(ias, sizeof(GF_ELEMENT));
        for (j=0; j<ias; j++)
            T[i][j] = A[ctoo_r[ncolA-ias+i]][ctoo_c[ncolA-ias+j]];
    }
    // Start to diagonalize
    clock_t start_p, stop_p;
    start_p = clock();
    long long ops1=0;
    GF_ELEMENT quotient;
    nonzeros = 0;
    for (i=0; i<ncolA-ias; i++) {
        for (j=i+1; j<ncolA; j++) {
            if (A[ctoo_r[i]][ctoo_c[i]] == 0 && get_loglevel() == TRACE)
                printf("The diagonal element after re-ordering is nonzero.\n");
            // process the item on (j, i)
            if (A[ctoo_r[j]][ctoo_c[i]] != 0) {
                quotient = galois_divide(A[ctoo_r[j]][ctoo_c[i]], A[ctoo_r[i]][ctoo_c[i]]);
                ops1 += 1;
                // multiply-and-add the corresponding part in the inactive part
                if (j < ncolA-ias) {
                    // eliminating nonzeros in the first ncolA-ias rows (but below diagonal)
                    galois_multiply_add_region(U[j], U[i], quotient, ias);
                } else {
                    // eliminating nonzeros in the last ias rows
                    galois_multiply_add_region(T[j-(ncolA-ias)], U[i], quotient, ias);
                }
                ops1 += ias;    // This part of matrix in processing is lower triangular in part, so operations only needed in the back half (i.e., inactiavted part)
                // simultaneously do the same thing on right matrix B
                galois_multiply_add_region(B[ctoo_r[j]], B[ctoo_r[i]], quotient, ncolB);
                ops1 += ncolB;
                A[ctoo_r[j]][ctoo_c[i]] = 0;            // eliminate the item
                nonzeros += 1;
            }
        }
    }
    stop_p = clock();
    if (get_loglevel() == TRACE) {
        printf("Diagonalize active part (%lld nonzero elements) after inactivating %d took %.6f seconds, cost %lld operations.\n", nonzeros, ias, ((double) (stop_p-start_p))/CLOCKS_PER_SEC, ops1);
    }
    operations += ops1;

    // Copy U and T back, free U but keep T as it is still needed for further pivoting
    for (i=0; i<ncolA-ias; i++) {
        for (j=0; j<ias; j++)
            A[ctoo_r[i]][ctoo_c[ncolA-ias+j]] = U[i][j];
        free(U[i]);
    }
    free(U);
    for (i=0; i<ias; i++) {
        for (j=0; j<ias; j++)
            A[ctoo_r[ncolA-ias+i]][ctoo_c[ncolA-ias+j]] = T[i][j];
    }

    // Second round of pivoting
    // Zlatev pivoting on (ias x ias) matrix
    ssList *row_pivots_2nd = malloc(sizeof(ssList));
    ssList *col_pivots_2nd = malloc(sizeof(ssList));
    row_pivots_2nd->ssFirst = row_pivots_2nd->ssLast = NULL;
    col_pivots_2nd->ssFirst = col_pivots_2nd->ssLast = NULL;
    if (get_loglevel() == TRACE) {
        printf("Zlatev pivoting...\n");
        start_pivoting = clock();
    }
    int ias_2nd = zlatev_pivoting(ias, ias, T, row_pivots_2nd, col_pivots_2nd);
    if (get_loglevel() == TRACE) {
        stop_pivoting = clock();
        pivoting_time = ((double) (stop_pivoting - start_pivoting)) / CLOCKS_PER_SEC;
        printf("Zlatev pivoting consumed time: %.6f seconds.\n", pivoting_time);
    }
    // Update row/col mappings
    int *partial_r = (int *) calloc(ias, sizeof(int));  // partial arrays are allocated to avoid corrupting the original ctoo_ arrays
    int *partial_c = (int *) calloc(ias, sizeof(int));
    Subscript *spt_2nd_r = row_pivots_2nd->ssFirst;
    Subscript *spt_2nd_c = col_pivots_2nd->ssFirst;
    for (i=0; i<ias; i++) {
        partial_r[i] = ctoo_r[ncolA-ias+spt_2nd_r->index];
        partial_c[i] = ctoo_c[ncolA-ias+spt_2nd_c->index];
        spt_2nd_r = spt_2nd_r->next;
        spt_2nd_c = spt_2nd_c->next;
    }
    memcpy(&(ctoo_r[ncolA-ias]), partial_r, sizeof(int)*ias);
    memcpy(&(ctoo_c[ncolA-ias]), partial_c, sizeof(int)*ias);
    free_subscriptList(row_pivots_2nd);
    free_subscriptList(col_pivots_2nd);
    free(partial_r);
    free(partial_c);

    /* Perform forward substitution on the ias x ias dense inactivated matrix. */
    start_p = clock();
    /*
    if (get_loglevel() == TRACE) {
        nonzeros = 0;
        total = 0;
        for (i=0; i<ias; i++) {
            for (j=0; j<ias; j++) {
                total++;
                if (T[i][j] != 0)
                    nonzeros++;
            }
        }
        printf("Current percentage of nonzeros of the sub-matrix: %.2f%%\n", (double) nonzeros/total*100);
        printf("Perform forward substitution on the square sub-matrix of size: %d x %d\n", ias, ias);
    }
    */
    // Make a copy of the corresponding msg matrices of T before performing forward substitution. 
    GF_ELEMENT **msg_submatrix = calloc(ias, sizeof(GF_ELEMENT*));
    for (i=0; i<ias; i++){
        msg_submatrix[i] = calloc(ncolB, sizeof(GF_ELEMENT));
        for (j=0; j<ias; j++)
            T[i][j] = A[ctoo_r[ncolA-ias+i]][ctoo_c[ncolA-ias+j]];  // reuse allocated memory, but data needs refresh because ctoo_r/ctoo_c were updated
        memcpy(msg_submatrix[i], B[ctoo_r[ncolA-ias+i]], ncolB*sizeof(GF_ELEMENT));
    }

    long long ops = forward_substitute(ias, ias, ncolB, T, msg_submatrix);
    operations += ops;
    // Save the processed inactivated part back to A
    for (i=0; i<ias; i++) {
        for (j=0; j<ias; j++)
            A[ctoo_r[ncolA-ias+i]][ctoo_c[ncolA-ias+j]] = T[i][j];
        memcpy(B[ctoo_r[ncolA-ias+i]], msg_submatrix[i], ncolB*sizeof(GF_ELEMENT));
    }

    //free T, msg_submatrix
    for (i=0; i<ias; i++) {
        free(T[i]);
        free(msg_submatrix[i]);
    }
    free(T);
    free(msg_submatrix);
    stop_p = clock();
    if (get_loglevel() == TRACE) {
        printf("Forward substitution on the bottom-right part took %.6f seconds, cost %lld operations\n", ((double) (stop_p-start_p))/CLOCKS_PER_SEC, ops);
    }

    /*
    if (get_loglevel() == TRACE) {
        int missing_pivots = 0;
        for (i=0; i<ncolA; i++) {
            if (A[ctoo_r[i]][ctoo_c[i]] == 0)
                missing_pivots += 1;
        }
        printf("There are %d pivots missing after forward substitution\n", missing_pivots);
    }
    */
    return operations;
}


/*****************************************************************************
 *     Zlatev pivoting
 *
 * A special kind of Markowitz pivoting in which pivots are selected from 3
 * candidates who have the smallest nonzeros.
 ******************************************************************************/
static int zlatev_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots)
{
    // Nonzeros of each rows and cols
    int *row_counts = (int *) calloc(nrow, sizeof(int));
    int *col_counts = (int *) calloc(ncolA, sizeof(int));

    int i, j, k;
    int max_row1s = 0;          // Max number of nonzeros in a row
    int max_col1s = 0;          // Max number of nonzeros in a column
    for (i=0; i<nrow; i++) {
        for (j=0; j<ncolA; j++) {
            if (A[i][j] != 0) {
                row_counts[i] += 1;
                if (row_counts[i] > max_row1s)
                    max_row1s = row_counts[i];

                col_counts[j] += 1;
                if (col_counts[j] > max_col1s)
                    max_col1s = col_counts[j];
            }
        }
    }
    /*************************************************************************
     * A list of double lists. Elements of each double list contain rows/cols
     * indices. Rows/cols of same double-linked list have the same number of
     * nonzero elements.
     **************************************************************************/
    ssList **RowID_lists = (ssList **) malloc(sizeof(ssList*) * (max_row1s+1));
    ssList **ColID_lists = (ssList **) malloc(sizeof(ssList*) * (max_col1s+1));
    for (i=0; i<max_row1s+1; i++) {
        RowID_lists[i] = (ssList *) malloc(sizeof(ssList));
        RowID_lists[i]->ssFirst = RowID_lists[i]->ssLast = NULL;
    }
    for (i=0; i<max_col1s+1; i++) {
        ColID_lists[i] = (ssList *) malloc(sizeof(ssList));
        ColID_lists[i]->ssFirst = ColID_lists[i]->ssLast = NULL;
    }
    // Travel through the matrix and populate double-linked lists
    int allzero_rows = 0;
    for (i=0; i<nrow; i++) {
        int row_count = row_counts[i];
        Subscript *rowID = malloc(sizeof(Subscript));
        rowID->index = i;
        rowID->nonzeros = row_count;
        rowID->prev = rowID->next = NULL;
        insertSubAtBeginning(RowID_lists[row_count], rowID);
        if (row_count == 0)
            allzero_rows += 1;
    }
    int allzero_cols = 0;
    for (i=0; i<ncolA; i++) {
        int col_count = col_counts[i];
        Subscript *colID = malloc(sizeof(Subscript));
        colID->index = i;
        colID->nonzeros = col_count;
        colID->prev = colID->next = NULL;
        insertSubAtBeginning(ColID_lists[col_count], colID);
        if (col_count == 0)
            allzero_cols += 1;
    }

    int nonallzero_rows = nrow - allzero_rows;
    int nonallzero_cols = ncolA - allzero_cols;
    if (get_loglevel() == TRACE)
        printf("There are %d all-zero rows and %d all-zero cols in the matrix.\n", allzero_rows, allzero_cols);
    // Markowitz pivoting using row_counts, col_counts, RowID_lists, ColID_lists
    int pivots_found = 0;
    int removed_cols = 0;
    int toallzero_cols = 0;
    int removed_rows = 0;
    int toallzero_rows = 0;
    int singletons = 0;
    while (pivots_found != ncolA) {
        //printf("%d rows removed, %d reduced to all-zero.\n", removed_rows, toallzero_rows);
        //printf("%d cols removed, %d reduced to all-zero.\n", removed_cols, toallzero_cols);

        // 计算Markowitz count
        int potential_r  = -1;
        int potential_c  = -1;
        int potential_mc = -1;              // potential Markowitz count
        int current_mc;
        int mc_minipos;             // the minimum possible Markowitz count in each row test
        int direct_found = 0;
        // search for one pivot
        Subscript *rows_inchecking, *cols_inchecking;
        int row_id, col_id;

        // Search rows
        int searched_rows = 0;
        for (i=1; i<=max_row1s; i++) {
            // 先按行的非零元素多少搜索
            mc_minipos = (i - 1) * (i - 1);
            if (RowID_lists[i]->ssFirst == NULL)
                continue;

            rows_inchecking = RowID_lists[i]->ssFirst;
            while (rows_inchecking != NULL && (searched_rows < ZLATEVS)) {
                searched_rows += 1;
                row_id = rows_inchecking->index;
                // 再按列的非零元素多少做test
                for (j=1; j<=max_col1s; j++) {
                    cols_inchecking = ColID_lists[j]->ssFirst;
                    while (cols_inchecking != NULL) {
                        col_id = cols_inchecking->index;
                        if (A[row_id][col_id] != 0) {
                            // we have found an entry in (row_id)-th row
                            current_mc = (i-1) * (j-1);
                            if (current_mc == 0) {

                                potential_r = row_id;
                                potential_c = col_id;
                                if (i==1)
                                    singletons += 1;
                                //  printf("a singleton row is found.\n");
                                goto found;
                            } else if (potential_mc == -1 || current_mc < potential_mc) {
                                potential_r = row_id;
                                potential_c = col_id;
                                potential_mc = current_mc;
                            }
                        }
                        cols_inchecking = cols_inchecking->next;
                    }
                }
                rows_inchecking = rows_inchecking->next;
            }
        }

        if (potential_r == -1 || potential_c == -1) {
            if (get_loglevel() == TRACE) {
                printf("%d rows reduced to all-zero.\n", toallzero_rows);
                printf("%d cols reduced to all-zero.\n", toallzero_cols);
                printf("(partial success) %d/%d pivots were found out of %d rows.\n", pivots_found, ncolA, nrow);
            }
            // There are row/col being reduced to all-zero, take them
            // as pivot anyway because we have no other choices
            Subscript *sub_ptt;
            sub_ptt = ColID_lists[0]->ssFirst;
            int zerocols = 0;
            while(sub_ptt != NULL) {
                zerocols += 1;
                Subscript *newCpivot0 = malloc(sizeof(Subscript));
                newCpivot0->index = sub_ptt->index;
                newCpivot0->nonzeros = sub_ptt->nonzeros;
                newCpivot0->prev = newCpivot0->next = NULL;
                insertSubAtEnd(ColPivots, newCpivot0);
                sub_ptt = sub_ptt->next;
            }

            sub_ptt = RowID_lists[0]->ssFirst;
            //for (i=0; i<(toallzero_rows+allzero_cols); i++)
            for (i=0; i<zerocols; i++) {
                Subscript *newRpivot0 = malloc(sizeof(Subscript));
                newRpivot0->index = sub_ptt->index;
                newRpivot0->nonzeros = sub_ptt->nonzeros;
                newRpivot0->prev = newRpivot0->next = NULL;
                insertSubAtEnd(RowPivots, newRpivot0);
                sub_ptt = sub_ptt->next;
            }
            if (get_loglevel() == TRACE)
                printf("There are %d/%d singleton rows were found as pivots.\n", singletons, pivots_found);
            return ncolA;
        }

found:
        // Found a pivot, save it
        pivots_found += 1;
        int p_r = potential_r;
        int p_c = potential_c;
        // Save coordiate (p_r, p_c)
        Subscript *newRpivot = malloc(sizeof(Subscript));
        newRpivot->index = p_r;
        newRpivot->nonzeros = row_counts[p_r];
        newRpivot->prev = newRpivot->next = NULL;
        insertSubAtEnd(RowPivots, newRpivot);
        Subscript *newCpivot = malloc(sizeof(Subscript));
        newCpivot->index = p_c;
        newCpivot->nonzeros = col_counts[p_c];
        newCpivot->prev = newCpivot->next = NULL;
        insertSubAtEnd(ColPivots, newCpivot);

        // Update row_counts[], col_counts[], RowID_lists, ColID_lists
        // 1, check nonzero elements of the row: potential_r, update numbers of nonzero elements of the correspondings columns(col_counts[] and ColID_lists).
        int nzs;
        Subscript *ss_pt;
        Subscript *ss_pt_next;
        // Note: some row/col have been eliminated, so we need to traverse Subscript when updating
        // 1, Update columns
        for (i=1; i<=max_col1s; i++) {
            ss_pt = ColID_lists[i]->ssFirst;
            while (ss_pt != NULL) {
                if (A[p_r][ss_pt->index] != 0) {
                    // 该subscript对象将要被更新处理，因此要记下当前遍历的位置
                    ss_pt_next = ss_pt->next;
                    if (ss_pt->index == p_c) {
                        removeSubscript(ColID_lists[i], ss_pt);
                        removed_cols += 1;
                        col_counts[ss_pt->index] -= 1;
                        free(ss_pt);
                    } else {
                        removeSubscript(ColID_lists[i], ss_pt);
                        ss_pt->nonzeros -= 1;
                        insertSubAtBeginning(ColID_lists[ss_pt->nonzeros], ss_pt);
                        if (ss_pt->nonzeros == 0)
                            toallzero_cols += 1;
                        col_counts[ss_pt->index] -= 1;
                    }
                    ss_pt = ss_pt_next;
                } else {
                    ss_pt = ss_pt->next;
                }
            }
        }
        // 2, Update rows
        for (j=1; j<=max_row1s; j++) {
            ss_pt = RowID_lists[j]->ssFirst;
            while (ss_pt != NULL) {
                if (A[ss_pt->index][p_c] != 0) {
                    // 该subscript对象将要被更新处理，因此要记下当前遍历的位置
                    ss_pt_next = ss_pt->next;
                    if (ss_pt->index == p_r) {
                        removeSubscript(RowID_lists[j], ss_pt);
                        removed_rows += 1;
                        row_counts[ss_pt->index] -= 1;
                        free(ss_pt);
                    } else {
                        removeSubscript(RowID_lists[j], ss_pt);
                        ss_pt->nonzeros -= 1;
                        if (ss_pt->nonzeros == 0)
                            toallzero_rows += 1;
                        insertSubAtBeginning(RowID_lists[ss_pt->nonzeros], ss_pt);
                        row_counts[ss_pt->index] -= 1;
                    }
                    ss_pt = ss_pt_next;
                } else {
                    ss_pt = ss_pt->next;
                }
            }
        }

    }

    free(row_counts);
    free(col_counts);
    for (i=0; i<max_row1s+1; i++)
        free_subscriptList(RowID_lists[i]);
    for (i=0; i<max_col1s+1; i++)
        free_subscriptList(ColID_lists[i]);
    if (get_loglevel() == TRACE) {
        printf("%d rows reduced to all-zero.\n", toallzero_rows);
        printf("%d cols reduced to all-zero.\n", toallzero_cols);
        printf("(full success) %d/%d pivots were found out of %d rows.\n", pivots_found, ncolA, nrow);
        printf("There are %d/%d singleton rows were found as pivots.\n", singletons, pivots_found);
    }
    return pivots_found;
}

/*********************************************************************************************************
 *      Inactivation pivoting
 * Use progressive inactivation to do pivoting. Only select pivots from singleton rows of the "residual"
 * matrix during the pivoting. If no singleton rows can be found, inactivate some columns until singletons
 * can be found:
 *  1) inactivate some columns and only perform pivoting on the rest of the "active" sub-matrix
 *  2) if singleton row cannot be found in the middle of pivoting, declare more inactive columns
 *  3) given the structure (heavier columns are in the back), declare inactive columns from the back
 *********************************************************************************************************/
static int inactivation_pivoting(int nrow, int ncolA, GF_ELEMENT **A, ssList *RowPivots, ssList *ColPivots)
{
    if (get_loglevel() == TRACE)
        printf("Pivoting matrix of size %d x %d via inactivation.\n", nrow, ncolA);

    int i, j, k;
    // 对矩阵中初始非零元素进行计数
    int *row_counts = (int *) calloc(nrow, sizeof(int));
    int *col_counts = (int *) calloc(ncolA, sizeof(int));

    int max_row1s = 0;                  // 记录初始矩阵里行中非零元素数目的最大值
    int max_col1s = 0;                  // 记录初始矩阵里列中非零元素数目的最大值
    for (i=0; i<nrow; i++) {
        for (j=0; j<ncolA; j++) {
            if (A[i][j] != 0) {
                row_counts[i] += 1;
                if (row_counts[i] > max_row1s)
                    max_row1s = row_counts[i];

                col_counts[j] += 1;
                if (col_counts[j] > max_col1s)
                    max_col1s = col_counts[j];
            }
        }
    }
    // 创建指针数组，每个指针指向一个双向链表，同一个链表中的列具有相同数目的非零元素
    // 该链表用来保存active的列的标号，按非零元素个数排列是为了方便inactivate含非零元素最多的列
    ssList **ColID_lists = (ssList **) malloc(sizeof(ssList*) * (max_col1s+1));
    for (i=0; i<max_col1s+1; i++) {
        ColID_lists[i] = (ssList *) malloc(sizeof(ssList));
        ColID_lists[i]->ssFirst = ColID_lists[i]->ssLast = NULL;
    }
    for (i=0; i<ncolA; i++) {
        int col_count = col_counts[i];
        Subscript *colID = malloc(sizeof(Subscript));
        colID->index = i;
        colID->nonzeros = col_count;
        colID->prev = colID->next = NULL;
        insertSubAtBeginning(ColID_lists[col_count], colID);
    }

    // Declare an array to record three possible "states" of each column
    // 0 - active
    // 1 - inactivated
    // 2 - removed because an entry of it was chosen as the pivot
    // Initially, all columns are active; in the end of the process, columns are
    // either inactivated or chosen as pivots in the active part
    uint8_t *col_state = (uint8_t *) calloc(ncolA, sizeof(uint8_t));

    int inactivated = 0;        // number of inactivated columns
    int active = ncolA;         // number of active columns

    int p_r;                    // used to store row index of the chosen pivot
    int p_c;                    // used to store col index of the chosen pivot
    int singleton_r_found;
    int ia;
    int selected_pivots = 0;
    while (active != 0) {
        p_r = -1;
        p_c = -1;
        singleton_r_found = 0;
        for (i=0; i<nrow; i++) {
            if (row_counts[i] == 1) {
                singleton_r_found = 1;
                //printf("singleton row is found at row %d.\n", i);
                p_r = i;
                break;
            }
        }

        if (singleton_r_found == 1) {
            // a singleton row is found, store the pivot
            for (j=0; j<ncolA; j++) {
                if ((col_state[j]==0) && (A[p_r][j]!=0)) {
                    p_c = j;
                    break;
                }
            }
            if (p_c == -1) {
                printf("error: failed to find the nonzero element in the singlton row.\n");
                exit(1);
            }
            // save the pivot
            // 保存该pivot的坐标
            Subscript *newRpivot = malloc(sizeof(Subscript));
            newRpivot->index = p_r;
            newRpivot->nonzeros = 1;
            newRpivot->prev = newRpivot->next = NULL;
            insertSubAtEnd(RowPivots, newRpivot);
            Subscript *newCpivot = malloc(sizeof(Subscript));
            newCpivot->index = p_c;
            newCpivot->nonzeros = col_counts[p_c];
            newCpivot->prev = newCpivot->next = NULL;
            insertSubAtEnd(ColPivots, newCpivot);

            // 更新row_counts
            row_counts[p_r] = -1;           // use -1 to indicate the row has an elelemnt was chosen as pivot
            for (i=0; i<nrow; i++) {
                if ( (row_counts[i] != -1) && (A[i][p_c] != 0) ) {
                    row_counts[i] -= 1;
                }
            }
            // 更新ColID_list
            Subscript *ss_pt = ColID_lists[col_counts[p_c]]->ssFirst;
            while (ss_pt != NULL) {
                if (ss_pt->index == p_c)
                    break;
                ss_pt = ss_pt->next;
            }
            if (ss_pt->index == p_c) {
                removeSubscript(ColID_lists[col_counts[p_c]], ss_pt);
                free(ss_pt);
                ss_pt = NULL;
            }

            col_state[p_c] = 2;
            active -= 1;
            selected_pivots += 1;
            //printf("Inactivated: %d Pivots Found: %d\n", inactivated, selected_pivots);
        } else {
            // no singleton row can be found, declare one column with the most nonzeros as inactive
            // Note: other algorithms may be used to choose a column to inactivate
            for (i=max_col1s; i>=0; i--) {
                Subscript *ss_pt = ColID_lists[i]->ssFirst;
                //if ((ss_pt != NULL) && (col_state[ss_pt->index] == 0))
                if (ss_pt != NULL) {
                    col_state[ss_pt->index] = 1;
                    inactivated += 1;
                    active -= 1;
                    for (k=0; k<nrow; k++) {
                        if ((A[k][ss_pt->index] != 0) && (row_counts[k] != -1))
                            row_counts[k] -= 1;
                    }

                    removeSubscript(ColID_lists[i], ss_pt);
                    free(ss_pt);
                    ss_pt = NULL;
                    break;
                }
            }
            selected_pivots += 1;
            //printf("Inactivated: %d Pivots Found: %d\n", inactivated, selected_pivots);
        }

    }
    // assign pivots for the dense part (any ordering is fine)
    for (i=0; i<ncolA; i++) {
        if (col_state[i] == 1) {
            // find an arbitray row that has not been selected
            int j_candidate;
            for (j=0; j<nrow; j++) {
                if (row_counts[j] != -1) {
                    j_candidate = j;
                    if (A[j][i] != 0)
                        break;
                }
            }
            j = j_candidate;

            // 保存该pivot的坐标
            Subscript *newRpivot = malloc(sizeof(Subscript));
            newRpivot->index = j;
            newRpivot->nonzeros = 0;
            newRpivot->prev = newRpivot->next = NULL;
            insertSubAtEnd(RowPivots, newRpivot);
            Subscript *newCpivot = malloc(sizeof(Subscript));
            newCpivot->index = i;
            newCpivot->nonzeros = 0;
            newCpivot->prev = newCpivot->next = NULL;
            insertSubAtEnd(ColPivots, newCpivot);
            row_counts[j] = -1;
            col_state[i] = 2;
            //printf("pivot at (%d, %d) is chosen for the dense part.\n", j, i);
        }
    }

    // free up memories
    free(col_state);

    for (i=0; i<max_col1s+1; i++)
        free_subscriptList(ColID_lists[i]);
    free(ColID_lists);

    free(row_counts);
    free(col_counts);
    return inactivated;
}

// insert a subscript object at the beginning of the list
static void insertSubAtBeginning(ssList *sub_list, Subscript *sub)
{
    if (sub_list->ssFirst == NULL) {
        sub_list->ssFirst = sub_list->ssLast = sub;
        sub->prev = NULL;
        sub->next = NULL;
    } else {
        sub_list->ssFirst->prev = sub;
        sub->next = sub_list->ssFirst;
        sub_list->ssFirst = sub;
        sub->prev = NULL;
    }
}

// insert a subscript object at the end of the list
static void insertSubAtEnd(ssList *sub_list, Subscript *sub)
{
    if (sub_list->ssFirst == NULL && sub_list->ssLast == NULL) {
        sub_list->ssFirst = sub_list->ssLast = sub;
        sub->prev = NULL;
        sub->next = NULL;
    } else {
        sub_list->ssLast->next = sub;
        sub->prev = sub_list->ssLast;
        sub->next = NULL;
        sub_list->ssLast = sub;
    }
}

// remove a subscript object from the list (no free() operation yet)
static void removeSubscript(ssList *sub_list, Subscript *sub)
{
    if (sub->prev == NULL && sub->next == NULL) {
        // only one left in the list
        sub_list->ssFirst = sub_list->ssLast = NULL;
    } else if (sub->prev == NULL && sub->next != NULL) {
        // the element is the head
        sub_list->ssFirst = sub->next;
        sub_list->ssFirst->prev = NULL;
    } else if (sub->prev != NULL && sub->next == NULL) {
        // the element is the tail
        sub_list->ssLast = sub->prev;
        sub->prev->next = NULL;
    } else {
        // the element is in the middle of the list
        sub->next->prev = sub->prev;
        sub->prev->next = sub->next;
    }
}

// free the subscript list
static void free_subscriptList(ssList *sub_list)
{
    Subscript *sub = sub_list->ssFirst;
    while (sub != NULL) {
        sub_list->ssFirst = sub->next;
        free(sub);
        sub = sub_list->ssFirst;
    }
    free(sub_list);
    sub_list = NULL;
}


