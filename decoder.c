#include <string.h>
#include "galois.h"
#include "bipartite.h"
#include "raptor.h"
static int apply_parity_check_matrix(struct dec_context *dc);
static void finish_recovering_BD(struct dec_context *dc);

static int apply_precode_matrix_and_pivot(struct dec_context *dc);
static void diag_decoding_matrix(struct dec_context *dc);

extern long long forward_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long long back_substitute(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B);
extern long pivot_matrix_oneround(int nrow, int ncolA, int ncolB, GF_ELEMENT **A, GF_ELEMENT **B, int **ctoo_r, int **ctoo_c, int *inactives);

// create decoding context
struct dec_context *create_decoder_context(struct enc_context *sc)
{
    int i, j, k;

    struct dec_context *dc = calloc(1, sizeof(struct dec_context));
    dc->sc = sc;

    dc->finished     = 0;
    dc->received     = 0;
    dc->deprecode    = 0;
    dc->inactivated  = 0;
    dc->dof          = 0;

    int snum    = dc->sc->snum;
    int pktsize = dc->sc->psize;
    int numpp   = dc->sc->snum + dc->sc->cnum;

    // allocate memory for A and B of Ax=B
    dc->coefficient = calloc(numpp, sizeof(GF_ELEMENT*));
    dc->message     = calloc(numpp, sizeof(GF_ELEMENT*));
    for (i=0; i<numpp; i++) {
        dc->coefficient[i] = calloc(numpp, sizeof(GF_ELEMENT));
        dc->message[i]     = calloc(pktsize, sizeof(GF_ELEMENT));
    }

    dc->ctoo_r = malloc(sizeof(int) * numpp);
    dc->ctoo_c = malloc(sizeof(int) * numpp);

    dc->operations   = 0;
    return dc;
}


void process_LT_packet(struct dec_context *dc, struct LT_packet *pkt)
{
    // dc->received += 1;

    int i, j, k;
    GF_ELEMENT quotient;

    int pktsize = dc->sc->psize;
    int numpp   = dc->sc->snum + dc->sc->cnum;


    if (dc->deprecode == 0) {
        // Before precode's check matrix was applied.
        // Simply save the coding vector in A of Ax=B at the corresponding row
        for (i=0; i<pkt->deg; i++) {
            dc->coefficient[dc->received][pkt->sid[i]] = 1;
        }
        memcpy(dc->message[dc->received], pkt->syms, pktsize*sizeof(GF_ELEMENT));
        dc->received += 1;
        if (dc->received == dc->sc->snum) {
            // Apply precode matrix to the last cnum rows of A, and inactivation pivoting
            int missing_DoF = apply_precode_matrix_and_pivot(dc);
            dc->dof = numpp - missing_DoF;
            dc->deprecode = 1;
        }
    } else {
        // Precode matrix has been applied, so the decoding matrix has been pivoted and re-ordered.
        // Process the received EV against the decoding matrix according to the new pivot order 
        dc->received += 1;
        GF_ELEMENT *ces = calloc(numpp, sizeof(GF_ELEMENT));
        for (i=0; i<pkt->deg; i++) {
                ces[pkt->sid[i]] = 1;
        }

        int pivot = -1;
        for (int m=0; m<numpp; m++) {
            if (ces[dc->ctoo_c[m]] != 0) {
                if (dc->coefficient[dc->ctoo_r[m]][dc->ctoo_c[m]] != 0) {
                    // mask the encoding vector and message over the JMB decoding matrix
                    GF_ELEMENT quotient = galois_divide(ces[dc->ctoo_c[m]], dc->coefficient[dc->ctoo_r[m]][dc->ctoo_c[m]]);
                    dc->operations += 1;
                    for (j=m; j<numpp; j++) {
                        ces[dc->ctoo_c[j]] = galois_add(ces[dc->ctoo_c[j]], galois_multiply(dc->coefficient[dc->ctoo_r[m]][dc->ctoo_c[j]], quotient));
                    }
                    galois_multiply_add_region(pkt->syms, dc->message[dc->ctoo_r[m]], quotient, pktsize);
                    dc->operations += pktsize;
                } else {
                    pivot = m;
                    break;
                }
            }
        }
        if (pivot != -1) {
            memcpy(dc->coefficient[dc->ctoo_r[pivot]], ces, numpp*sizeof(GF_ELEMENT));
            memcpy(dc->message[dc->ctoo_r[pivot]], pkt->syms,  pktsize*sizeof(GF_ELEMENT));
            dc->dof += 1;
        }
        free(ces);
        ces = NULL;
    }


    if (dc->dof == numpp) {
        diag_decoding_matrix(dc);
    }

}

// Apply the parity-check matrix to the decoding matrix and perform inactivation pivoting
static int apply_precode_matrix_and_pivot(struct dec_context *dc)
{
    int i, j, k;
    int num_of_new_DoF = 0;

    int pktsize = dc->sc->psize;
    int numpp = dc->sc->snum + dc->sc->cnum;

    // 1, Copy parity-check vectors to the nonzero rows of the decoding matrix
    // Apply precoding matrix
    for (i=0; i<dc->sc->cnum; i++) {
        dc->coefficient[dc->sc->snum+i][dc->sc->snum+i] = 1;
        NBR_node *variable_node = dc->sc->graph->l_nbrs_of_r[i]->first;  //ldpc_graph->nbrs_of_right[i];
        while (variable_node != NULL) {
            // go through source packet nodes connecting to this check node
            int src_pktid = variable_node->data;                        //variable_node->nb_index;
            dc->coefficient[dc->sc->snum+i][src_pktid] = variable_node->ce;
            variable_node = variable_node->next;
        }
    }

    // 2, Pivot decoding matrix, and return reorder sequence as arguments
    dc->operations += pivot_matrix_oneround(numpp, numpp, pktsize, dc->coefficient, dc->message, &dc->ctoo_r, &dc->ctoo_c, &(dc->inactivated));

    /* Count available innovative rows */
    int missing_DoF = 0;
    for (i=0; i<numpp; i++) {
        if (dc->coefficient[dc->ctoo_r[i]][dc->ctoo_c[i]] == 0)
            missing_DoF++;
    }
    return missing_DoF;
}


static void diag_decoding_matrix(struct dec_context *dc)
{
    int pktsize = dc->sc->psize;
    int numpp = dc->sc->snum + dc->sc->cnum;
    int i, j, k;
    GF_ELEMENT quotient;
    long long bs_ops = 0;
    // Backard substitution from right-most col to the left
    for (j=numpp-1; j>=0; j--) {
        for (i=0; i<j; i++) {
            if (dc->coefficient[dc->ctoo_r[i]][dc->ctoo_c[j]] != 0) {
                quotient = galois_divide(dc->coefficient[dc->ctoo_r[i]][dc->ctoo_c[j]], dc->coefficient[dc->ctoo_r[j]][dc->ctoo_c[j]]);
                galois_multiply_add_region(dc->message[dc->ctoo_r[i]], dc->message[dc->ctoo_r[j]], quotient, pktsize);
                dc->coefficient[dc->ctoo_r[i]][dc->ctoo_c[j]] = 0;
                bs_ops += 1 + pktsize;
            }
        }
    }
    // Convert all diagonal element to 1
    dc->pp = calloc(numpp, sizeof(GF_ELEMENT*));
    for (i=0; i<numpp; i++) {
        if (dc->coefficient[dc->ctoo_r[i]][dc->ctoo_c[i]] != 1)
            galois_multiply_region(dc->message[dc->ctoo_r[i]], galois_divide(1, dc->coefficient[dc->ctoo_r[i]][dc->ctoo_c[i]]), pktsize);
        dc->coefficient[dc->ctoo_r[i]][dc->ctoo_c[i]] = 1;
        // recover decoded data packets
        int pktid = dc->ctoo_c[i];
        dc->pp[pktid] = calloc(pktsize, sizeof(GF_ELEMENT));
        memcpy(dc->pp[pktid], dc->message[dc->ctoo_r[i]], pktsize*sizeof(GF_ELEMENT));
    }
    dc->operations += bs_ops;
    dc->finished = 1;
}

void free_decoder_context(struct dec_context *dc)
{
    if (dc == NULL)
        return;
    if (dc->coefficient != NULL) {
        for (int i=dc->sc->snum+dc->sc->cnum-1; i>=0; i--) {
            if (dc->coefficient[i] != NULL)
                free(dc->coefficient[i]);
        }
        free(dc->coefficient);
    }
    if (dc->message != NULL) {
        for (int i=dc->sc->snum+dc->sc->cnum-1; i>=0; i--) {
            if (dc->message[i] != NULL)
                free(dc->message[i]);
        }
        free(dc->message);
    }
    if (dc->pp != NULL) {
        for (int i=dc->sc->snum+dc->sc->cnum-1; i>=0; i--) {
            if (dc->pp[i] != NULL)
                free(dc->pp[i]);
        }
        free(dc->pp);
    }
    if (dc->ctoo_r != NULL)
        free(dc->ctoo_r);
    if (dc->ctoo_c != NULL)
        free(dc->ctoo_c);
    free(dc);
    dc = NULL;
    return;
}
