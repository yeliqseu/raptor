#include <string.h>
#include "raptor.h"
#include "galois.h"
#include "bipartite.h"

static double dist[40];    // degree distribution for LT encoding

static int compare_int(const void *elem1, const void *elem2);
static void precoding(struct enc_context *sc);
static int draw_random_degree(void);
static void get_random_unique_numbers(int ids[], int n, int ub);

// Create an encoder packet from a buf of data
struct enc_context *create_encoder_context(GF_ELEMENT *buf, int snum, int pktsize)
{
    struct enc_context *sc = malloc(sizeof(struct enc_context));
    sc->snum  = snum;
    sc->psize = pktsize;
    sc->cnum  = 13;            // use fixed number for now
    sc->count = 0;

    construct_GF();

    // load data
    int i;
    sc->pp = calloc(sc->snum+sc->cnum, sizeof(GF_ELEMENT*));
    for (i=0; i<sc->snum; i++) {
        sc->pp[i] = calloc(sc->psize, sizeof(GF_ELEMENT));
        memcpy(sc->pp[i], buf, sc->psize*sizeof(GF_ELEMENT));
        buf += sc->psize;
    }

    // allocate parity-check packet space
    for (i=0; i<sc->cnum; i++)
        sc->pp[sc->snum+i] = calloc(sc->psize, sizeof(GF_ELEMENT));
    
    // precoding
    sc->graph = malloc(sizeof(BP_graph));
    sc->graph->binaryce = 0;                   // binary or non-binary precoding
    create_bipartite_graph(sc->graph, sc->snum, sc->cnum); 
    precoding(sc);

    // configure LT degree distribution
    // ref. R10 code
    dist[0]  = 0.00971;
    dist[1]  = 0.4580;
    dist[2]  = 0.2100;
    dist[3]  = 0.1130;
    dist[9]  = 0.1110;
    dist[10] = 0.0797;
    dist[39] = 0.0156;
    return sc;
}

// Perform systematic LDPC precoding 
static void precoding(struct enc_context *sc)
{
    int i, j;
    for (i=0; i<sc->cnum; i++) {
        // Encoding check packet according to the LDPC graph
        NBR_node *nb = sc->graph->l_nbrs_of_r[i]->first;
        while(nb != NULL) {
            int sid = nb->data;  // index of source packet
            // XOR information content
            galois_multiply_add_region(sc->pp[i+sc->snum], sc->pp[sid], nb->ce, sc->psize);
            // move to next possible neighbour node of current check
            nb = nb->next;
        }
    }
}

// Encode an LT packet from the intermediate packets
struct LT_packet *encode_LT_packet(struct enc_context *sc)
{
    struct LT_packet *pkt = calloc(1, sizeof(struct LT_packet));
    pkt->id = sc->count;
    int deg = draw_random_degree();
    if (deg > sc->snum+sc->cnum) {
        deg = sc->snum+sc->cnum;        // for small number of packets, the largest degree might be too large
    }
    pkt->deg = deg;
    pkt->sid = calloc(pkt->deg, sizeof(int));
    // draw source packet id
    get_random_unique_numbers(pkt->sid, pkt->deg, sc->snum+sc->cnum);
    // combine packets 
    pkt->syms = calloc(sc->psize, sizeof(GF_ELEMENT));
    for (int i=0; i<pkt->deg; i++) {
        galois_multiply_add_region(pkt->syms, sc->pp[pkt->sid[i]], 1, sc->psize);
    }
    sc->count += 1;
    return pkt;
}

struct LT_packet *duplicate_LT_packet(struct LT_packet *spkt, struct enc_context *sc)
{
    struct LT_packet *pkt = calloc(1, sizeof(struct LT_packet));
    pkt->id = spkt->id;
    pkt->deg = spkt->deg;
    pkt->sid = calloc(pkt->deg, sizeof(int));
    pkt->syms = calloc(sc->psize, sizeof(GF_ELEMENT));
    memcpy(pkt->sid, spkt->sid, sizeof(int)*pkt->deg);
    memcpy(pkt->syms, spkt->syms, sizeof(GF_ELEMENT)*sc->psize);
    return pkt;
}

int free_LT_packet(struct LT_packet *pkt)
{
    free(pkt->syms);
    free(pkt->sid);
    free(pkt);
    return (0);
}


void free_encoder_context(struct enc_context *sc)
{
    if (sc == NULL)
        return;
    int i;
    if (sc->pp != NULL) {
        for (i=sc->snum+sc->cnum-1; i>=0; i--) {
            if (sc->pp[i] != NULL) {
                free(sc->pp[i]);
                sc->pp[i] = NULL;
            }
        }
        free(sc->pp);
    }
    if (sc->graph != NULL)
        free_bipartite_graph(sc->graph);
    free(sc);
    sc = NULL;
    return;
}

static int draw_random_degree(void)
{
    int ds = 1;
    int dm = 40;
    int Precision = 100000;
    int degree = 1;
    double degree_int[dm-ds+1];
    for (int i=0; i<dm-ds+1; i++)
        degree_int[i] = dist[i] * Precision;

    int choice = rand() % Precision;
    double lower_bound = 0;
    double upper_bound = lower_bound + degree_int[0];
    for (int j=0; j<dm-ds+1; j++) {
        if (choice >= lower_bound && choice < upper_bound) {
            degree = ds + j;
            break;
        } else {
            lower_bound += degree_int[j];
            upper_bound = lower_bound + degree_int[j+1];
        }
    }
    return degree;
}

// generate a number of n<ub unique random numbers within the range of [0, ub-1]
// using Fisher-Yates shuffle method
static void get_random_unique_numbers(int ids[], int n, int ub)
{
	int init_array[ub];
	int i, j;
	for (i=0; i<ub; i++)
		init_array[i] = i;

	// randomly shuffle the init_array
	for (i=ub-1; i>=1; i--) {
		int rd = rand() % (i+1);
		//int rd = gsl_rng_uniform_int(r, i+1);
		int tmp = init_array[rd];
		init_array[rd] = init_array[i];
		init_array[i] = tmp;
	}

    // sort the obtained unique random numbers so that coding coefficients corresponding
    // to packets are stored in the ascending order (to simplify decoder implementation)
    qsort(init_array, n, sizeof(int), compare_int);
    memcpy(ids, init_array, n*sizeof(int));
	//for (j=0; j<n; j++)
	//	ids[j] = init_array[j];
}

static int compare_int(const void *elem1, const void *elem2)
{
    int a = * ((int *) elem1);
    int b = * ((int *) elem2);
    if (a < b)
        return -1;
    if (a > b)
        return 1;
    return 0;
}
