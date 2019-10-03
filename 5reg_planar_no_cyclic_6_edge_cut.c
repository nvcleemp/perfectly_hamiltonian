#include <stdlib.h>
#include <stdio.h>
#include "planegraphs_base.h"
#include "planegraphs_input.h"
#include "planegraphs_output.h"

typedef struct __cycle_data {
    boolean *current_cycle;
} CYCLE_DATA;

boolean continue_cycle(PLANE_GRAPH *graph, int last, int remaining, int first, CYCLE_DATA *data) {
    if(remaining==0){
        //no more vertices left
        PG_EDGE *e, *e_last;
        e = e_last = graph->first_edge[last];
        do {
            if(e->end == first){
                return TRUE;
            }
            e = e->next;
        } while (e != e_last);
        return FALSE;
    } else {
        //try to extend the cycle
        PG_EDGE *e, *e_last;
        e = e_last = graph->first_edge[last];
        do {
            if(!data->current_cycle[e->end]){
                data->current_cycle[e->end] = TRUE;
                if(continue_cycle(graph, e->end, remaining-1, first, data)){
                    return TRUE;
                }
                data->current_cycle[e->end] = FALSE;
            }
            e = e->next;
        } while (e != e_last);
    }

    return FALSE;
}

/**
 * Finds a cycle of length target_length containing start_vertex
 * @param graph
 * @param start_vertex
 * @param target_length
 * @param data
 * @return TRUE if there is a cycle of length target_length containing start_vertex, and FALSE otherwise
 */
boolean start_cycle(PLANE_GRAPH *graph, int start_vertex, int target_length, CYCLE_DATA *data){
    //mark the start vertex as being in the cycle
    data->current_cycle[start_vertex] = TRUE;

    //go through neighbours of start vertex and start the search for a cycle
    //we don't need to try the last neighbour, since the cycle needs to return to this vertex
    PG_EDGE *e, *e_last;
    e = graph->first_edge[start_vertex];
    e_last = e->prev;
    do {
        if(!data->current_cycle[e->end]){
            data->current_cycle[e->end] = TRUE;
            //search for cycle containing the edge e
            if(continue_cycle(graph, e->end, target_length - 2, start_vertex, data)){
                return TRUE;
            }
            data->current_cycle[e->end] = FALSE;
        }
        e = e->next;
    } while (e != e_last);

    data->current_cycle[start_vertex] = FALSE;

    return FALSE;
}


boolean has_six_cycle(PLANE_GRAPH *graph){
    CYCLE_DATA data;
    data.current_cycle = malloc(sizeof(boolean)*graph->nv);

    for (int i = 0; i < graph->nv; ++i) {
        if(graph->degree[i] > 3){
            //set current cycle to FALSE, except for vertices that already were start vertices
            for (int j = 0; j < graph->nv; ++j) {
                data.current_cycle[j] = (j < i) && (graph->degree[j] > 3);
            }
            if(start_cycle(graph, i, 6, &data)){
                free(data.current_cycle);
                return TRUE;
            }
        }
    }

    free(data.current_cycle);
    return FALSE;
}

int main(int argc, char *argv[]) {
    unsigned long graph_count = 0;
    unsigned long written_count = 0;

    DEFAULT_PG_INPUT_OPTIONS(options);
    PLANE_GRAPH *graph;
    while((graph = read_and_decode_planar_code(stdin, &options))){
        graph_count++;
        for (int i = 0; i < graph->nv; ++i) {
            if(graph->degree[i] != 5){
                fprintf(stderr, "Graph %lu is not a quintic graph -- exiting!\n", graph_count);
                exit(EXIT_FAILURE);
            }
        }

        PLANE_GRAPH *dual = get_dual_graph(graph);
        if(!has_six_cycle(dual)){
            written_count++;
            write_planar_code(graph, stdout);
        }

        free_plane_graph(dual);
        free_plane_graph(graph);
    }

    fprintf(stderr, "Read %lu graph%s.\n", graph_count, graph_count==1 ? "" : "s");
    fprintf(stderr, "Written %lu graph%s without cyclic-6-edge-cut.\n", written_count, written_count==1 ? "" : "s");
}