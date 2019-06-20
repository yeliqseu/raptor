/*----------------------- bipartite.c ----------------------
 *
 *  Bipartite graph functions for LDPC precoding.
 *  We use LDPC code specified in Raptor Codes Standard
 *  Implementation. Refer to:
 *    1, Sec 5.4.2.3 in RFC5053 "Raptor Forward Error
 *       Correction Scheme for Object Delivery" by Luby et. al.
 *    2, Sec 3.2.3 in Foundations and Trends in Comm. and
 *       Info. Theory "Raptor Codes" by Shokrollahi et. al.
 *  for detailed information.
 *
 *----------------------------------------------------------*/
#include "bipartite.h"
#include "galois.h"
#include <math.h>
#include <string.h>
static int is_prime(int number);
static int include_left_node(int l_index, int r_index, BP_graph *graph);
static void free_list(struct node_list *list);
static void clear_list(struct node_list *list);
static int exist_in_list(struct node_list *list, int data);
static int remove_from_list(struct node_list *list, int data);
static void append_to_list(struct node_list *list, struct node *nd);

// construct LDPC graph
int create_bipartite_graph(BP_graph *graph, int nleft, int nright)
{
    int LDPC_SYS = nleft;
    int S        = nright;
    if (S == 0)
        return 0;

    int i, j;
    graph->nleft  = nleft;
    graph->nright = nright;
    if ((graph->l_nbrs_of_r = calloc(nright, sizeof(NBR_nodes*))) == NULL)
        goto failure;
    for (i=0; i<nright; i++) {
        if ((graph->l_nbrs_of_r[i] = malloc(sizeof(NBR_nodes))) == NULL)
            goto failure;
        graph->l_nbrs_of_r[i]->first = graph->l_nbrs_of_r[i]->last = NULL;
    }

    if ((graph->r_nbrs_of_l = calloc(nleft, sizeof(NBR_nodes*))) == NULL)
        goto failure;
    for (i=0; i<nleft; i++) {
        if ((graph->r_nbrs_of_l[i] = malloc(sizeof(NBR_nodes))) == NULL)
            goto failure;
        graph->r_nbrs_of_l[i]->first = graph->r_nbrs_of_l[i]->last = NULL;
    }

    char *hdpc = getenv("SNC_PRECODE");
    if (hdpc != NULL && strcmp(hdpc, "HDPC") == 0) {
        // A reference bipartite graph, which is highly dense. This is only used 
        // when SNC_PRECODE env var is set to HDPC. The env var is ONLY for development
        // and testing use.
        for (i=0; i<S; i++) {
            for (j=0; j<LDPC_SYS; j++) {
                int included = 1;
                if (graph->binaryce == 1) {
                    if (rand() % 2 == 0)
                        included = 0;
                } else {
                    if (rand() % 256 == 0)
                        included = 0;
                }
                if (included) {
                    if (include_left_node(j, i, graph) < 0)
                        goto failure;
                }
            }
        }
        return 0;
    }

    // By default use circulant LDPC code which is used by Raptor code
    int a, b;
    int touching_edge = 0;
    for (i=0; i<ceil((double) LDPC_SYS/S); i++) {
        if (touching_edge == 1)
            break;

        // assign non-zero positions for the first column in each circulant matrix
        // each check node connects to exactly 3 left nodes
        // 1, P[0][i*S] = 1;
        if (include_left_node(i*S, 0, graph) < 0)
            goto failure;
        //2, P[a-1][i*S] = 1;
        a = (((i+1)+1)%S == 0) ? S : ((i+1)+1)%S;
        if (include_left_node(i*S, a-1, graph) < 0)
            goto failure;
        //3, P[b-1][i*S] = 1;
        b = ((2*(i+1)+1)%S == 0) ? S : (2*(i+1)+1)%S;
        if (include_left_node(i*S, b-1, graph) < 0)
            goto failure;

        // circulant part
        for (j=1; j<S; j++) {
            if (i*S + j >= LDPC_SYS) {
                touching_edge = 1;
                break;
            }
            // shift down the non-zero positions of previous columns in the circulant matrix
            //1, P[0+j][i*S+j] = 1;
            if (include_left_node(i*S+j, 0+j, graph) < 0)
                goto failure;
            //2, P[a-1][i*S+j] = 1;
            a = (((i+1)+1+j)%S == 0) ? S : ((i+1)+1+j)%S;
            if (include_left_node(i*S+j, a-1, graph) < 0)
                goto failure;
            //3, P[b-1][i*S+j] = 1;
            b = ((2*(i+1)+1+j)%S == 0) ? S : (2*(i+1)+1+j)%S;
            if (include_left_node(i*S+j, b-1, graph) < 0)
                goto failure;
        }
    }
    return 0;

failure:
    free_bipartite_graph(graph);
    return -1;
}

// include left node index in the LDPC graph
static int include_left_node(int l_index, int r_index, BP_graph *graph)
{
    // Skip if the two nodes are already neighbors
    // Note: a good ``bipartitin'' algorithm should not get into such 
    // trouble. This is included just in case we needed to test/experiment
    // different bipartite creating methods.
    if (exist_in_list(graph->l_nbrs_of_r[r_index], l_index))
        return 0;
    // Coding coefficient associated with the edge
    GF_ELEMENT ce;
    if (graph->binaryce == 1) {
        ce = 1;
    } else {
        ce = (GF_ELEMENT) (rand() % 255 + 1); // Value range: [1-255]
    }
    // Record neighbor of a right-side node
    NBR_node *nb = calloc(1, sizeof(NBR_node));
    if (nb == NULL)
        return -1;
    nb->data = l_index;
    nb->ce   = ce;
    nb->next = NULL;
    append_to_list(graph->l_nbrs_of_r[r_index], nb);

    // Also neighbour of the left-side node
    NBR_node *check_nb = calloc(1, sizeof(NBR_node));
    if (nb == NULL)
        return -1;
    check_nb->data = r_index;
    check_nb->ce   = ce;
    check_nb->next = NULL;
    append_to_list(graph->r_nbrs_of_l[l_index], check_nb);
    return 0;
}

void free_bipartite_graph(BP_graph *graph)
{
    if (graph == NULL)
        return;
    int i;
    if (graph->r_nbrs_of_l != NULL) {
        for (i=0; i<graph->nleft; i++)
            free_list(graph->r_nbrs_of_l[i]);
        free(graph->r_nbrs_of_l);
    }
    if (graph->l_nbrs_of_r != NULL) {
        for (i=0; i<graph->nright; i++)
            free_list(graph->l_nbrs_of_r[i]);
        free(graph->l_nbrs_of_r);
    }
    free(graph);
}


static void append_to_list(struct node_list *list, struct node *nd)
{
    if (list->first == NULL)
        list->first = list->last = nd;
    else {
        list->last->next = nd;
        list->last = nd;
    }
}

// Remove the first node whose data is equal to "data"
// Note: this function should only be used in applications
//       where nodes in the list have unique data
static int remove_from_list(struct node_list *list, int data)
{
    struct node *prev = NULL;
    struct node *curr = list->first;
    while (curr != NULL) {
        if (curr->data == data) {
            // shorten list
            if (curr == list->first && curr == list->last) {            // list contains only one node
                list->first = list->last = NULL;
            } else if (curr == list->first && curr != list->last) {     // head node is to be removed
                list->first = curr->next;
            } else if (curr != list->first && curr == list->last) {     // tail node is to be removed
                list->last = prev;
                list->last->next = NULL;
            } else {
                prev->next = curr->next;
            }

            free(curr);
            curr = NULL;
            return 0;
        }
        prev = curr;
        curr = curr->next;
    }
    return -1;
}

static int exist_in_list(struct node_list *list, int data)
{
    struct node *p = list->first;
    while (p != NULL) {
        if (p->data == data)
            return 1;
        p = p->next;
    }
    return 0;
}

// clear nodes in a list, but keep the list structure alive
static void clear_list(struct node_list *list)
{
    struct node *nd = list->first;
    struct node *ndNext;
    while (nd != NULL) {
        ndNext = nd->next;
        free(nd);
        nd = ndNext;
    }
    list->first = list->last = NULL;
}

// Free a list, which include clear nodes in a list and free
// the list structure in the end.
static void free_list(struct node_list *list)
{
    if (list != NULL) {
        clear_list(list);
        free(list);
    }
}

