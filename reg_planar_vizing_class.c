#include <stdlib.h>
#include <stdio.h>
#include <bitset.h>
#include <getopt.h>
#include <unionfind.h>
#include "planegraphs_base.h"
#include "planegraphs_input.h"
#include "planegraphs_output.h"

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

boolean complete_colouring(PLANE_GRAPH *graph){
    int next_vertex = 1; //0 is guaranteed to be fully coloured
    while(next_vertex < graph->nv && colours_at_vertex[next_vertex]==fully_coloured){
        next_vertex++;
    }

    if(next_vertex==graph->nv){
        //graph is fully coloured
        return TRUE;
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

                if(complete_colouring(graph)){
                    return TRUE;
                }

                //revert colour
                e->colour = e->inverse->colour = -1;
                NEIGHBOUR_ALONG_COLOUR(e->start, i) = -1;
                NEIGHBOUR_ALONG_COLOUR(e->end, i) = -1;
                BS_REMOVE(colours_at_vertex[e->start], i);
                BS_REMOVE(colours_at_vertex[e->end], i);
            }
        }
    }
    return FALSE;
}

boolean is_class_1(PLANE_GRAPH *graph){
    PG_EDGE *e;

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

    boolean result = complete_colouring(graph);

    free(neighbour_along_colour);
    free(colours_at_vertex);

    return result;
}

//====================== USAGE =======================

void help(char *name) {
    fprintf(stderr, "This program determines the Vizing class for k-regular plane graphs.\n\n");
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options] k\n\n", name);
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "    -v, --verbose\n");
    fprintf(stderr, "       Print the result for each graph to standard error.\n");
    fprintf(stderr, "    -1, --one\n");
    fprintf(stderr, "       Filter the graphs that are Vizing class 1.\n");
    fprintf(stderr, "    -2, --two\n");
    fprintf(stderr, "       Filter the graphs that are Vizing class 2.\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options] k\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    boolean verbose = FALSE;
    boolean one = FALSE;
    boolean two = FALSE;

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
            {"verbose", no_argument, NULL, 'v'},
            {"one", no_argument, NULL, '1'},
            {"two", no_argument, NULL, '2'},
            {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'v':
                verbose = TRUE;
                break;
            case '1':
                one = TRUE;
                break;
            case '2':
                two = TRUE;
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

    if(one && two){
        fprintf(stderr, "Class 1 and class 2 graphs cannot be exported at the same time -- exiting!\n");
        return EXIT_FAILURE;
    }

    if(!(one || two)){
        verbose = TRUE; //otherwise there will be no output.
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
    int written_count = 0;
    while((graph = read_and_decode_planar_code(stdin, &options))){
        graph_count++;
        for (int i = 0; i < graph->nv; ++i) {
            if(graph->degree[i] != k){
                fprintf(stderr, "Not a k-regular graph -- exiting!\n");
                exit(EXIT_FAILURE);
            }
        }

        if(is_class_1(graph)){
            if(verbose) {
                fprintf(stderr, "Graph %d is class 1.\n", graph_count);
            }
            if(one) {
                written_count++;
                write_planar_code(graph, stdout, written_count == 1);
            }
        } else {
            if (verbose) {
                fprintf(stderr, "Graph %d is class 2.\n", graph_count);
            }
            if(two) {
                written_count++;
                write_planar_code(graph, stdout, written_count == 1);
            }
        }

        free_plane_graph(graph);
    }

    fprintf(stderr, "Read %d graph%s.\n\n", graph_count, graph_count==1 ? "" : "s");
    if(one || two) {
        fprintf(stderr, "Written %d graph%s that have Vizing class %d.\n\n", written_count,
                written_count == 1 ? "" : "s", one ? 1 : 2);
    }

}