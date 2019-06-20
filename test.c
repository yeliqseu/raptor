#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "raptor.h"

char usage[] = "usage: ./raptorTest snum nhop epsils\n\
                       snum    - number of source packets\n\
                       nhop    - number of hops\n\
                       epsils - erasure probability of each hop (one number if all hops are the same)\n";

int main(int argc, char *argv[])
{
    if (argc <4 || (argc !=4 && argc != atoi(argv[2])+3)) {
        printf("%s\n", usage);
        exit(1);
    }
    int i, j;
    // read network configurations
    int snum = atoi(argv[1]);
    int psize = 200;
    int nhop = atoi(argv[2]);
    double *pe = calloc(nhop, sizeof(double));
    for (i=0; i<nhop; i++) {
        if (argc == 4) {
            pe[i] = atof(argv[3]);
        } else {
            pe[i] = atof(argv[3+i]);
        }
        printf("pe[%d]: %.7f ", i, pe[i]);
    }
    printf("\n");
    // construct encoder, decoder, and relay buffers
    struct timeval tv;
    gettimeofday(&tv, NULL);
    srand(tv.tv_sec * 1000 + tv.tv_usec / 1000); // seed use microsec
    unsigned char *buf = malloc(snum*psize);
    int rnd=open("/dev/urandom", O_RDONLY);
    read(rnd, buf, snum*psize);
    close(rnd);

    // create encoder context
    struct enc_context *sc = create_encoder_context(buf, snum, psize);
    // create buffers
    struct LT_buffer *rbuf = calloc(nhop-1, sizeof(struct LT_buffer));
    for (i=0; i<nhop-1; i++) {
        rbuf[i].pkt = NULL;
    }

    // create decoder context
    struct dec_context *dc = create_decoder_context(sc);

    int nuse = 0;
    struct LT_packet *pkt = NULL;
    while (!dc->finished) {
        nuse++;
        for (i=0; i<nhop; i++) {
            // Transmission of the hop
            if (i==0) {
                // first hop. generate a packet
                pkt = encode_LT_packet(sc);
            } else {
                // intermediate node forwarding if they have buffered packets
                if (rbuf[i-1].pkt != NULL) {
                    // pkt = duplicate_LT_packet(rbuf[i-1].pkt, sc);
                    //free_LT_packet(rbuf[i-1].pkt);   // simulating 0-length buffer
                    //rbuf[i-1].pkt = NULL;
                    pkt = rbuf[i-1].pkt;    // simulating 0-length buffer
                    rbuf[i-1].pkt = NULL;
                } else {
                    continue;  // no packet sent, so skip the receiving part
                }
            }
            // Receiving of the hop
            if (rand() % 10000000 >= pe[i]*10000000) {
                // successfully received
                if (i<nhop - 1) {
                    // intermediate node
                    rbuf[i].pkt = pkt;
                } else {
                    // destination
                    process_LT_packet(dc, pkt);
                    free_LT_packet(pkt);
                }
            } else {
                free_LT_packet(pkt);
                pkt = NULL;
            }
        }
    }
    for (i=0; i<snum; i++) {
        if (memcmp(sc->pp[i], dc->pp[i], sizeof(GF_ELEMENT)*psize) != 0)
            fprintf(stderr, "recovered packet %d is NOT identical to original.\n", i);
    }
    printf("\npktsize: %d nhop: %d snum: %d cnum: %d received: %d nuse: %d \n", psize, nhop, snum, sc->cnum, dc->received, nuse);
    free_decoder_context(dc);
    free(buf);
    free_encoder_context(sc);
    return 0;
}
