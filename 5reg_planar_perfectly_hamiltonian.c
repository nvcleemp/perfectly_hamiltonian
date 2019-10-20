#include <stdlib.h>
#include <stdio.h>
#include <bitset.h>
#include <getopt.h>
#include "planegraphs_base.h"
#include "planegraphs_input.h"

//Stores the neighbour of a vertex along a specific colour.
//Either the label of the neighbour is stored or -1 if this colour is not
//yet assigned at that vertex.
//Best to address this array using the macro below.
int *neighbour_along_colour;
#define NEIGHBOUR_ALONG_COLOUR(v, c) neighbour_along_colour[5*(v) + c]

bitset *colours_at_vertex;

bitset fully_coloured = BS_SINGLETON(5) - 1ULL;

boolean print_perfectly_hamiltonian_colouring = FALSE;

void print_colouring(PLANE_GRAPH *graph){
    for (int i = 0; i < graph->nv; ++i) {
        fprintf(stdout, "%2d: ", i);
        for (int j = 0; j < 5; ++j) {
            fprintf(stdout, " %2d", NEIGHBOUR_ALONG_COLOUR(i, j));
        }
        fprintf(stdout, "\n");
    }
}

int length_of_two_colour_cycle(int start, int colour1, int colour2){
    int next = NEIGHBOUR_ALONG_COLOUR(start, colour2);
    int length = 1;
    while(TRUE){//ugly but we clean it up later (or not)
        if(next==-1) return -1;
        next = NEIGHBOUR_ALONG_COLOUR(next, colour1);
        length++;
        if(next==-1) return -1;
        if(next==start) return length;
        next = NEIGHBOUR_ALONG_COLOUR(next, colour2);
        length++;
    }
}

boolean complete_colouring(PLANE_GRAPH *graph){
    int next_vertex = 1; //0 is guaranteed to be fully coloured
    while(next_vertex < graph->nv && colours_at_vertex[next_vertex]==fully_coloured){
        next_vertex++;
    }

    if(next_vertex==graph->nv){
        //graph is fully coloured and all cycles have already been checked
        if(print_perfectly_hamiltonian_colouring) print_colouring(graph);
        return TRUE;
    } else {
        //choose an edge to be coloured
        PG_EDGE *e = graph->first_edge[next_vertex];
        while(e->index!=-1) e = e->next; //will stop because next_vertex is not fully coloured

        //determine possible colours
        for (int i = 0; i < 5; ++i) {
            if(!BS_CONTAINS(BS_UNION(colours_at_vertex[e->start], colours_at_vertex[e->end]),i)){
                //assign colour
                e->index = e->inverse->index = i;
                NEIGHBOUR_ALONG_COLOUR(e->start, i) = e->end;
                NEIGHBOUR_ALONG_COLOUR(e->end, i) = e->start;
                BS_ADD(colours_at_vertex[e->start], i);
                BS_ADD(colours_at_vertex[e->end], i);

                boolean colours_still_OK = TRUE;
                for (int j = 0; j < 5; ++j) {
                    if(j!=i){
                        int length = length_of_two_colour_cycle(e->end, i, j);
                        if(length > 0 && length < graph->nv){
                            colours_still_OK = FALSE;
                            break;
                        }
                    }
                }

                if(colours_still_OK && complete_colouring(graph)){
                    return TRUE;
                }

                //revert colour
                e->index = e->inverse->index = -1;
                NEIGHBOUR_ALONG_COLOUR(e->start, i) = -1;
                NEIGHBOUR_ALONG_COLOUR(e->end, i) = -1;
                BS_REMOVE(colours_at_vertex[e->start], i);
                BS_REMOVE(colours_at_vertex[e->end], i);
            }
        }
        return FALSE;
    }
}

boolean is_perfectly_hamiltonian(PLANE_GRAPH *graph){
    PG_EDGE *e;

    neighbour_along_colour = (int *)malloc(5*graph->nv* sizeof(int));
    colours_at_vertex = (bitset *)malloc(graph->nv* sizeof(bitset));
    for (int i = 0; i < graph->nv; ++i) {
        e = graph->first_edge[i];
        for (int j = 0; j < 5; ++j) {
            NEIGHBOUR_ALONG_COLOUR(i, j) = -1;
            e->index = -1;
            e = e->next;
        }
        colours_at_vertex[i] = BS_EMPTY_SET;
    }

    //colour edges at first vertex
    e = graph->first_edge[0];
    for (unsigned int i = 0; i < 5; ++i) {
        NEIGHBOUR_ALONG_COLOUR(0, i) = e->end;
        NEIGHBOUR_ALONG_COLOUR(e->end, i) = 0;
        BS_ADD(colours_at_vertex[0], i);
        BS_ADD(colours_at_vertex[e->end], i);
        e->index = e->inverse->index = i;
        e = e->next;
    }

    boolean result = complete_colouring(graph);

    free(neighbour_along_colour);
    free(colours_at_vertex);
    return result;
}


//====================== USAGE =======================

void help(char *name) {
    fprintf(stderr, "This program finds perfectly hamiltonian 5-regular plane graphs.\n\n");
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -p, --print\n");
    fprintf(stderr, "       Print the colours for the perfectly hamiltonian graphs.\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options]\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
            {"print", no_argument, NULL, 'p'},
            {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "hp", long_options, &option_index)) != -1) {
        switch (c) {
            case 'p':
                print_perfectly_hamiltonian_colouring = TRUE;
                break;
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

    /*=========== main loop ===========*/

    DEFAULT_PG_INPUT_OPTIONS(options);
    PLANE_GRAPH *graph;
    int graph_count = 0;
    while((graph = read_and_decode_planar_code(stdin, &options))){
        graph_count++;
        for (int i = 0; i < graph->nv; ++i) {
            if(graph->degree[i] != 5){
                fprintf(stderr, "Not a quintic graph -- exiting!\n");
                exit(EXIT_FAILURE);
            }
        }

        if(is_perfectly_hamiltonian(graph)){
            fprintf(stderr, "Graph %d is perfectly hamiltonian.\n", graph_count);
        } else {
            fprintf(stderr, "Graph %d is not perfectly hamiltonian.\n", graph_count);
        }
        free_plane_graph(graph);
    }
}