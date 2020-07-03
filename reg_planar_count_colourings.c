#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include <bitset.h>

#include <planegraphs_base.h>
#include <planegraphs_input.h>

//the degree of the graph: will be set to either 3, 4, or 5
int k;

//Stores the neighbour of a vertex along a specific colour.
//Either the label of the neighbour is stored or -1 if this colour is not
//yet assigned at that vertex.
//Best to address this array using the macro below.
int *neighbour_along_colour;
#define NEIGHBOUR_ALONG_COLOUR(v, c) neighbour_along_colour[k*(v) + c]

bitset *colours_at_vertex;

bitset fully_coloured;

int current_count = 0;

void handle_colouring(PLANE_GRAPH *graph){
    current_count++;
}

void complete_colouring(PLANE_GRAPH *graph){
    int next_vertex = 1; //0 is guaranteed to be fully coloured
    while(next_vertex < graph->nv && colours_at_vertex[next_vertex]==fully_coloured){
        next_vertex++;
    }

    if(next_vertex==graph->nv){
        //graph is fully coloured
        handle_colouring(graph);
    } else {
        //choose an edge to be coloured
        PG_EDGE *e = graph->first_edge[next_vertex];
        while(e->colour!=-1) e = e->next; //will stop since next_vertex is not fully coloured

        //determine possible colours
        for (int i = 0; i < k; ++i) {
            if(!BS_CONTAINS(BS_UNION(colours_at_vertex[e->start], colours_at_vertex[e->end]),i)){
                //assign colour
                e->colour = e->inverse->colour = i;
                NEIGHBOUR_ALONG_COLOUR(e->start, i) = e->end;
                NEIGHBOUR_ALONG_COLOUR(e->end, i) = e->start;
                BS_ADD(colours_at_vertex[e->start], i);
                BS_ADD(colours_at_vertex[e->end], i);

                complete_colouring(graph);

                //revert colour
                e->colour = e->inverse->colour = -1;
                NEIGHBOUR_ALONG_COLOUR(e->start, i) = -1;
                NEIGHBOUR_ALONG_COLOUR(e->end, i) = -1;
                BS_REMOVE(colours_at_vertex[e->start], i);
                BS_REMOVE(colours_at_vertex[e->end], i);
            }
        }
    }
}

void all_colourings(PLANE_GRAPH *graph){
    PG_EDGE *e;

    current_count = 0;

    neighbour_along_colour = (int *)malloc(k*graph->nv* sizeof(int));
    colours_at_vertex = (bitset *)malloc(graph->nv* sizeof(bitset));
    for (int i = 0; i < graph->nv; ++i) {
        e = graph->first_edge[i];
        for (int j = 0; j < k; ++j) {
            NEIGHBOUR_ALONG_COLOUR(i, j) = -1;
            e->colour = -1;
            e = e->next;
        }
        colours_at_vertex[i] = BS_EMPTY_SET;
    }

    //colour edges at first vertex
    e = graph->first_edge[0];
    for (unsigned int i = 0; i < k; ++i) {
        NEIGHBOUR_ALONG_COLOUR(0, i) = e->end;
        NEIGHBOUR_ALONG_COLOUR(e->end, i) = 0;
        BS_ADD(colours_at_vertex[0], i);
        BS_ADD(colours_at_vertex[e->end], i);
        e->colour = e->inverse->colour = i;
        e = e->next;
    }

    complete_colouring(graph);

    free(neighbour_along_colour);
    free(colours_at_vertex);
}

//====================== USAGE =======================

void help(char *name) {
    fprintf(stderr, "This program counts the number of 1-factorisations in regular plane graphs.\n\n");
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options] k\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
            {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'h':
                help(name);
                return EXIT_SUCCESS;
            case '?':
                usage(name);
                return EXIT_FAILURE;
            default:
                fprintf(stderr, "Illegal option %c.\n", c);
                usage(name);
                return EXIT_FAILURE;
        }
    }

    if (optind == argc) {
        usage(name);
        return EXIT_FAILURE;
    }

    k = atoi(argv[optind]);

    if(k < 3 || k > 5){
        fprintf(stderr, "k should be 3, 4, or 5 -- exiting!\n");
        exit(EXIT_FAILURE);
    }

    fully_coloured = BS_SINGLETON(k) - 1ULL;

    DEFAULT_PG_INPUT_OPTIONS(options);
    PLANE_GRAPH *graph;
    int graph_count = 0;
    while((graph = read_and_decode_planar_code(stdin, &options))){
        graph_count++;
        for (int i = 0; i < graph->nv; ++i) {
            if(graph->degree[i] != k){
                fprintf(stderr, "Not a k-regular graph -- exiting!\n");
                exit(EXIT_FAILURE);
            }
        }

        all_colourings(graph);

        fprintf(stderr, "Graph %d has %d edge-colourings.\n", graph_count, current_count);

        free_plane_graph(graph);
    }

    fprintf(stderr, "Read %d graph%s.\n\n", graph_count, graph_count==1 ? "" : "s");

}