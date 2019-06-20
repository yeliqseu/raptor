#ifndef GALOIS
#define GALOIS
typedef unsigned char GF_ELEMENT;
#endif

struct LT_packet {
    int         id;         // packet id
    int         deg;        // packet degree
    int         *sid;       // id of source packets in the packet 
    GF_ELEMENT  *syms;      //symbols of coded packet
};


struct enc_context {
    int                       snum;     // number of source packet
    int                       cnum;     // number of check packet
    int                       psize;    // packet size in bytes
    struct  bipartite_graph  *graph;    // bipartite graph of precoding
    GF_ELEMENT              **pp;       // pointers to precoded source packets
    int                       count;    // number of generated coded pacekts
};


// Raptor decoder context
struct dec_context
{
    struct enc_context *sc;

    int finished;               // an indicator tracking the finish of decoding
    int received;               // total received LT packet
    int deprecode;              // apply precode or not
    int inactivated;            // total number of inactivated packets
    int dof;                    // number of linearly independent rows

    // decoding matrix
    GF_ELEMENT **coefficient;
    GF_ELEMENT **message;

    // the following two mappings are to record pivoting processings
    int *ctoo_r;                // record the mapping from current row index to the original row id
    int *ctoo_c;                // record the mapping from current col index to the original row id

    GF_ELEMENT **pp;            // decoded packets
    /*performance index*/
    long long operations;       // record the number of computations used
};

// A buffer of size 0/1?
struct LT_buffer
{
    struct LT_packet *pkt;
};

// Encoder
struct enc_context *create_encoder_context(GF_ELEMENT *buf, int snum, int pktsize);
void free_encoder_context(struct enc_context *sc);
struct LT_packet *encode_LT_packet(struct enc_context *sc);
int free_LT_packet(struct LT_packet *pkt);

// Decoder
struct dec_context *create_decoder_context(struct enc_context *sc);
void process_LT_packet(struct dec_context *dc, struct LT_packet *pkt);
void free_decoder_context(struct dec_context *dc);
