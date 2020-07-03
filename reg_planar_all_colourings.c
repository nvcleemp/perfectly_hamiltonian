#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>

#include <bitset.h>
#include <unionfind.h>

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

//gives for each possible pair count the number of colourings obtaining this count
int pair_count_count[11];

int factorisation_count = 0;

typedef struct __1factorisation_list_element FACTORISATION_LIST_ELEMENT;

struct __1factorisation_list_element {
    bitset *factorisation;
    int index;
    int perfect_pair_count;

    FACTORISATION_LIST_ELEMENT *smaller;
    FACTORISATION_LIST_ELEMENT *larger;
};

FACTORISATION_LIST_ELEMENT *factorisation_tree;
FACTORISATION_LIST_ELEMENT **factorisation_list;

FACTORISATION_LIST_ELEMENT *new_factorisation_list_element(){
    FACTORISATION_LIST_ELEMENT *el = (FACTORISATION_LIST_ELEMENT *)malloc(sizeof(FACTORISATION_LIST_ELEMENT));
    if(el==NULL){
        fprintf(stderr, "Insufficient memory -- exiting!\n");
        exit(EXIT_FAILURE);
    }
    el->smaller = el->larger = NULL;
    el->factorisation = malloc(sizeof(bitset)*k);
    if(el->factorisation==NULL){
        fprintf(stderr, "Insufficient memory -- exiting!\n");
        exit(EXIT_FAILURE);
    }
    return el;
}

void free_factorisation_list_element(FACTORISATION_LIST_ELEMENT *el){
    if(el->smaller!=NULL) free_factorisation_list_element(el->smaller);
    if(el->larger!=NULL) free_factorisation_list_element(el->larger);
    free(el->factorisation);
    free(el);
}

void number_factorisation_list(FACTORISATION_LIST_ELEMENT *current, int *counter){
    if(current==NULL) return;
    number_factorisation_list(current->smaller, counter);
    current->index = *counter;
    (*counter)++;
    number_factorisation_list(current->larger, counter);
}

void build_factorisation_list(FACTORISATION_LIST_ELEMENT *current){
    if(current==NULL) return;
    factorisation_list[current->index] = current;
    build_factorisation_list(current->smaller);
    build_factorisation_list(current->larger);
}

int count_factorisation_list(FACTORISATION_LIST_ELEMENT *current){
    if(current==NULL) return 0;
    return count_factorisation_list(current->smaller) + count_factorisation_list(current->larger) + 1;
}

PG_EDGE *edges[BS_SET_CAPACITY];

void label_edges(PLANE_GRAPH *graph){
    if(graph->ne/2 > BS_SET_CAPACITY){
        fprintf(stderr, "Graph has too many edges -- exiting!\n");
        exit(EXIT_FAILURE);
    }

    RESETMARKS(graph);
    int edge_counter = 0;

    PG_EDGE *e, *elast;

    //write the original vertices: only adjacent to center of edges
    for(int i=0; i<graph->nv; i++){
        e = elast = graph->first_edge[i];
        do {
            if(!ISMARKED(graph, e)){
                e->index = e->inverse->index = edge_counter;
                edges[edge_counter] = e;
                edge_counter++;
                MARK(graph, e);
                MARK(graph, e->inverse);
            }
            e = e->next;
        } while (e != elast);
    }
}

void sort_vector_in_place(bitset *vector) {
    for(int i=1; i<k; i++){
        bitset temp = vector[i];
        int j=i-1;
        while(j>=0 && temp<vector[j]){
            vector[j+1] = vector[j];
            j--;
        }
        vector[j+1] = temp;
    }
}

void translate_to_1factor_vector(PLANE_GRAPH *graph, bitset *vector){
    for (int i = 0; i < k; ++i) {
        vector[i] = BS_EMPTY_SET;
        for (int j = 0; j < graph->ne / 2; ++j) {
            if(edges[j]->colour==i){
                BS_ADD(vector[i], j);
            }
        }
    }
    sort_vector_in_place(vector);
}

void apply_from_1factor_vector(PLANE_GRAPH *graph, bitset *vector){
    for (int i = 0; i < k; ++i) {
        for (int j = 0; j < graph->ne / 2; ++j) {
            if(BS_CONTAINS(vector[i], j)){
                edges[j]->colour=i;
                edges[j]->inverse->colour=i;
                NEIGHBOUR_ALONG_COLOUR(edges[j]->start, i) = edges[j]->end;
                NEIGHBOUR_ALONG_COLOUR(edges[j]->end, i) = edges[j]->start;
            }
        }
    }
}

void perform_kempe_switch(bitset component, bitset *old_factorization, bitset *new_factorization){
    for (int i = 0; i < k; ++i) {
        if(BS_IS_NOT_EMPTY(BS_INTERSECTION(component, old_factorization[i]))){
            new_factorization[i] = BS_SYMMETRIC_DIFFERENCE(old_factorization[i], component);
        } else {
            new_factorization[i] = old_factorization[i];
        }
    }
    sort_vector_in_place(new_factorization);
}

int get_components_for_1factor_pair(bitset f1, bitset f2, bitset *components){
    int component_count = 0;
    bitset uncovered_edges = f1;
    while(BS_IS_NOT_EMPTY(uncovered_edges)){
        bitset current_component = BS_EMPTY_SET;
        int start = __builtin_ffsll(uncovered_edges) - 1;
        PG_EDGE *start_edge = edges[start];
        PG_EDGE *current_edge = edges[start];
        do {
            BS_ADD(current_component, current_edge->index);
            current_edge = current_edge->inverse;
            while (!BS_CONTAINS(f2, current_edge->index)) current_edge = current_edge->next;
            BS_ADD(current_component, current_edge->index);
            current_edge = current_edge->inverse;
            while (!BS_CONTAINS(f1, current_edge->index)) current_edge = current_edge->next;
        } while (start_edge != current_edge);
        BS_SAFE_REMOVE_ALL(uncovered_edges, current_component);
        if(components!=NULL){
            components[component_count] = current_component;
        }
        component_count++;
    }
    return component_count;
}

void write_coloured_graph(FILE *f, PLANE_GRAPH *graph, int pair_count){
    static boolean first = TRUE;
    int order = graph->nv + graph->ne/2 + k + 1 + pair_count + 1;
    if(order > 253){
        fprintf(stderr, "Graphs of that size are currently not supported -- exiting!\n");
        return;
    }
    if(first){
        fprintf(f, ">>multi_code<<");
        first = FALSE;
    }

    PG_EDGE *e, *elast;

    //write the number of vertices
    fputc(order, f);

    //write the original vertices: only adjacent to center of edges
    int edge_center_start = graph->nv + 1;
    for(int i=0; i<graph->nv; i++){
        e = elast = graph->first_edge[i];
        do {
            fputc(edge_center_start + e->index, f);
            e = e->next;
        } while (e != elast);
        fputc(0, f);
    }
    //write centers of edges: only remaining neighbours are colours
    int colour_start = graph->nv + graph->ne/2 + 1;
    for(int i=0; i<graph->ne/2; i++){
        fputc(colour_start + edges[i]->colour, f);
        fputc(0, f);
    }
    for(int i=0; i<k; i++){
        fputc(colour_start + k, f);
        fputc(0, f);
    }
    for(int i=0; i<pair_count + 1; i++){
        fputc(colour_start + k + 1 + i, f);
    }
    fputc(0, f);
    for(int i=0; i<pair_count; i++){
        fputc(0, f);
    }
}

void print_colouring(FILE *f, PLANE_GRAPH *graph){
    for (int i = 0; i < graph->nv; ++i) {
        fprintf(f, "%2d: ", i);
        for (int j = 0; j < k; ++j) {
            fprintf(f, " %2d", NEIGHBOUR_ALONG_COLOUR(i, j));
        }
        fprintf(f, "\n");
    }
}

void print_colouring_from_vector(FILE *f, PLANE_GRAPH *graph, bitset *factorisation){
    int *temp = neighbour_along_colour;
    neighbour_along_colour = (int *)malloc(k*graph->nv* sizeof(int));
    apply_from_1factor_vector(graph, factorisation);
    print_colouring(f, graph);
    free(neighbour_along_colour);
    neighbour_along_colour = temp;
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

int compare_bitset_vector(bitset *vector1, bitset *vector2){
    for (int i = 0; i < k; ++i) {
        if(vector1[i] > vector2[i]) return 1;
        if(vector1[i] < vector2[i]) return -1;
    }
    return 0;
}

FACTORISATION_LIST_ELEMENT *find_1factorisation(bitset *vector){
    FACTORISATION_LIST_ELEMENT *current = factorisation_tree;
    int compare = compare_bitset_vector(vector, current->factorisation);
    while(compare){
        if(compare > 0){
            current = current->larger;
        } else {
            current = current->smaller;
        }
        compare = compare_bitset_vector(vector, current->factorisation);
    }
    return current;
}

void store_current_colouring(PLANE_GRAPH *graph, int perfect_pair_count){
    bitset vector[k];
    translate_to_1factor_vector(graph, vector);

    FACTORISATION_LIST_ELEMENT *current;
    if(factorisation_tree == NULL){
        current = factorisation_tree = new_factorisation_list_element();
    } else {
        current = factorisation_tree;
        while(TRUE){
            if(compare_bitset_vector(vector, current->factorisation) > 0){
                if(current->larger==NULL){
                    current->larger = new_factorisation_list_element();
                    current = current->larger;
                    break;
                }
                current = current->larger;
            } else {
                if(current->smaller==NULL){
                    current->smaller = new_factorisation_list_element();
                    current = current->smaller;
                    break;
                }
                current = current->smaller;
            }
        }
    }

    current->perfect_pair_count = perfect_pair_count;
    for (int i = 0; i < k; ++i) {
        current->factorisation[i] = vector[i];
    }
    factorisation_count++;
}

void handle_colouring(PLANE_GRAPH *graph){
    int good_pairs = 0;
    for (int i = 0; i < k-1; ++i) {
        for (int j = i+1; j < k; ++j) {
            if(length_of_two_colour_cycle(0, i, j) == graph->nv){
                good_pairs++;
            }
        }
    }
    pair_count_count[good_pairs]++;
    //write_coloured_graph(stdout, graph, good_pairs);
    store_current_colouring(graph, good_pairs);
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

                /*boolean colours_still_OK = TRUE;
                for (int j = 0; j < k; ++j) {
                    if(j!=i){
                        int length = length_of_two_colour_cycle(e->end, i, j);
                        if(length > 0 && length < graph->nv){
                            colours_still_OK = FALSE;
                            break;
                        }
                    }
                }*/

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
    fprintf(stderr, "This program counts 1-factorisations (or edge-colourings) in regular plane graphs.\n\n");
    fprintf(stderr, "Usage\n=====\n");
    fprintf(stderr, " %s [options]\n\n", name);
    fprintf(stderr, "Valid options\n=============\n");
    fprintf(stderr, "\nGeneral\n-------\n");
    fprintf(stderr, "    -c, --counts-per-graph\n");
    fprintf(stderr, "       Print the number of 1-factorisations per graph as well as\n");
    fprintf(stderr, "       the number of 1-factorisations for each number of perfect pairs.\n");
    fprintf(stderr, "    -K, --Kempe\n");
    fprintf(stderr, "       Compute the Kempe equivalence classes for each graph.\n");
    fprintf(stderr, "    -C, --classes-per-graph\n");
    fprintf(stderr, "       Give an overview of the Kempe equivalence classes per graph. This\n");
    fprintf(stderr, "       will print the number of classes, the size of each class and an\n");
    fprintf(stderr, "       overview of the number of perfect pair counts within this class.\n");
    fprintf(stderr, "       This will also set option -K.\n");
    fprintf(stderr, "    -h, --help\n");
    fprintf(stderr, "       Print this help and return.\n");
    fprintf(stderr, "\nKempe switches\n--------------\n");
    fprintf(stderr, "    -M, --maximum-from-single-switch\n");
    fprintf(stderr, "       Print for each graph the maximum number of different perfect pair\n");
    fprintf(stderr, "       counts that can be reached from a single 1-factorisation by performing\n");
    fprintf(stderr, "       only one switch.\n");
    fprintf(stderr, "    -F n, --from n\n");
    fprintf(stderr, "    -T k, --to k\n");
    fprintf(stderr, "       Print for each graph the 1-factorisations with n perfect pairs that\n");
    fprintf(stderr, "       can be transformed into a 1-factorisation with k perfect pairs and\n");
    fprintf(stderr, "       not in a 1-factorisation with a different number of perfect pairs.\n");
    fprintf(stderr, "       The option -T can be used several times, and it will increase the\n");
    fprintf(stderr, "       number of permitted reachable perfect pairs.\n");
    fprintf(stderr, "    -A, --at-least\n");
    fprintf(stderr, "       This option is used in combination with -F and -T, and the program will\n");
    fprintf(stderr, "       then also output 1-factorisations that do no reach all k's specified\n");
    fprintf(stderr, "       with the option -T.\n");
}

void usage(char *name) {
    fprintf(stderr, "Usage: %s [options] k\n", name);
    fprintf(stderr, "For more information type: %s -h \n\n", name);
}

int main(int argc, char *argv[]) {

    boolean counts_per_graph = FALSE;
    boolean classes_per_graph = FALSE;
    boolean maximum_from_single_switch = FALSE;
    int global_maximum_from_single_switch = 0; //will be at least 1
    int global_maximum_from_single_switch_per_count[11];
    for(int i = 0; i < 11; i++){
        global_maximum_from_single_switch_per_count[i] = -1;
    }

    boolean compute_equivalence_classes = FALSE;

    int from = -1;
    bitset target = BS_EMPTY_SET;
    boolean target_at_least = FALSE;

    /*=========== commandline parsing ===========*/

    int c;
    char *name = argv[0];
    static struct option long_options[] = {
            {"Kempe", no_argument, NULL, 'K'},
            {"kempe", no_argument, NULL, 'K'},
            {"counts-per-graph", no_argument, NULL, 'c'},
            {"classes-per-graph", no_argument, NULL, 'C'},
            {"maximum-from-single-switch", no_argument, NULL, 'M'},
            {"from", required_argument, NULL, 'F'},
            {"to", required_argument, NULL, 'T'},
            {"at-least", no_argument, NULL, 'A'},
            {"help", no_argument, NULL, 'h'}
    };
    int option_index = 0;

    while ((c = getopt_long(argc, argv, "KcCMF:T:Ah", long_options, &option_index)) != -1) {
        switch (c) {
            case 'K':
                compute_equivalence_classes = TRUE;
                break;
            case 'c':
                counts_per_graph = TRUE;
                break;
            case 'C':
                classes_per_graph = TRUE;
                compute_equivalence_classes = TRUE;
                break;
            case 'M':
                maximum_from_single_switch = TRUE;
                compute_equivalence_classes = TRUE;
                break;
            case 'F':
                from = atoi(optarg);
                compute_equivalence_classes = TRUE;
                break;
            case 'T':
                c = atoi(optarg);
                if(c >= 0 && c < BS_SET_CAPACITY){
                    BS_ADD(target, c);
                }
                compute_equivalence_classes = TRUE;
                break;
            case 'A':
                target_at_least = TRUE;
                compute_equivalence_classes = TRUE;
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

    if (optind == argc) {
        usage(name);
        return EXIT_FAILURE;
    }

    k = atoi(argv[optind]);

    if(k < 3 || k > 5){
        fprintf(stderr, "k should be 3, 4, or 5 -- exiting!\n");
        exit(EXIT_FAILURE);
    }

    boolean single_switch_attainable[k*(k-1)/2+1][k*(k-1)/2+1];
    for(int i = 0; i < k*(k-1)/2+1; i++){
        for(int j = 0; j < k*(k-1)/2+1; j++){
            single_switch_attainable[i][j] = FALSE;
        }
    }

    fully_coloured = BS_SINGLETON(k) - 1ULL;

    DEFAULT_PG_INPUT_OPTIONS(options);
    PLANE_GRAPH *graph;
    int graph_count = 0;
    int min_equivalence_class_count = -1;
    while((graph = read_and_decode_planar_code(stdin, &options))){
        factorisation_count = 0;
        graph_count++;
        for (int i = 0; i < graph->nv; ++i) {
            if(graph->degree[i] != k){
                fprintf(stderr, "Not a k-regular graph -- exiting!\n");
                exit(EXIT_FAILURE);
            }
        }
        if(graph->ne/2 > BS_SET_CAPACITY){
            fprintf(stderr, "Graph %d has too many edges -- skipping!\n", graph_count);
            continue;
        }
        label_edges(graph);

        for (int i = 0; i < k*(k-1)/2+1; ++i) {
            pair_count_count[i] = 0;
        }

        all_colourings(graph);

        if(counts_per_graph) {
            fprintf(stderr, "Graph %d:\n", graph_count);
            for (int i = 0; i < k * (k - 1) / 2 + 1; ++i) {
                fprintf(stderr, "  %2d: %6d\n", i, pair_count_count[i]);
            }
            fprintf(stderr, "%d 1-factorisation%s\n\n", factorisation_count, factorisation_count == 1 ? "" : "s");
        }

        if(compute_equivalence_classes){
            int factorisation_counter = 0;
            number_factorisation_list(factorisation_tree, &factorisation_counter);

            if(factorisation_counter != factorisation_count){
                fprintf(stderr, "Incorrect factorisation count -- exiting!\n");
                exit(EXIT_FAILURE);
            }

            factorisation_list = (FACTORISATION_LIST_ELEMENT **)malloc(sizeof(FACTORISATION_LIST_ELEMENT *)*factorisation_count);
            build_factorisation_list(factorisation_tree);

            int maximum_different_pairs_reachable = 0;
            int colouring_index = 0;
            //do Kempe switches and union-find
            UNIONFIND *uf_1factors = unionfind_prepare_new(factorisation_count);
            for (int i = 0; i < factorisation_count; ++i) {
                boolean reachable_from_current_with_single_switch[k*(k-1)/2+1];
                for(int j = 0; j < k*(k-1)/2+1; j++){
                    reachable_from_current_with_single_switch[j]=FALSE;
                }
                bitset reachable_bs = BS_EMPTY_SET;
                for (int j = 0; j < k - 1; ++j) {
                    for (int l = j + 1; l < k; ++l) {
                        bitset components[graph->nv/4];
                        int component_count = get_components_for_1factor_pair(
                                factorisation_list[i]->factorisation[j],
                                factorisation_list[i]->factorisation[l],
                                components);
                        for (int m = 0; m < component_count; ++m) {
                            bitset switched[k];
                            perform_kempe_switch(components[m], factorisation_list[i]->factorisation, switched);
                            FACTORISATION_LIST_ELEMENT *image = find_1factorisation(switched);
                            unionfind_union(uf_1factors, i, image->index);
                            single_switch_attainable[factorisation_list[i]->perfect_pair_count][image->perfect_pair_count] = TRUE;
                            reachable_from_current_with_single_switch[image->perfect_pair_count] = TRUE;
                            BS_ADD(reachable_bs, image->perfect_pair_count);
                        }
                    }
                }
                if(factorisation_list[i]->perfect_pair_count==from && (
                        reachable_bs==target ||
                        (target_at_least && BS_CONTAINS_ALL(reachable_bs, target)))){
                    fprintf(stderr, "From %d to ( ", from);
                    for(int j = 0; j < k*(k-1)/2+1; j++){
                        if(BS_CONTAINS(target, j)){
                            fprintf(stderr, "%d ", j);
                        }
                    }
                    fprintf(stderr, "):\n");
                    print_colouring_from_vector(stderr, graph, factorisation_list[i]->factorisation);
                }
                int maximum_different_pairs_reachable_from_current_colouring = 0;
                for(int j = 0; j < k*(k-1)/2+1; j++){
                    if(reachable_from_current_with_single_switch[j]){
                        maximum_different_pairs_reachable_from_current_colouring++;
                    }
                }
                if(global_maximum_from_single_switch_per_count[factorisation_list[i]->perfect_pair_count] < maximum_different_pairs_reachable_from_current_colouring){
                    global_maximum_from_single_switch_per_count[factorisation_list[i]->perfect_pair_count] = maximum_different_pairs_reachable_from_current_colouring;
                }
                if(maximum_different_pairs_reachable_from_current_colouring > maximum_different_pairs_reachable){
                    maximum_different_pairs_reachable = maximum_different_pairs_reachable_from_current_colouring;
                    colouring_index = i;
                }
            }

            if(global_maximum_from_single_switch < maximum_different_pairs_reachable){
                global_maximum_from_single_switch = maximum_different_pairs_reachable;
            }

            if(maximum_from_single_switch) {
                fprintf(stderr, "Maximum of different number of perfect pairs reachable from one colouring by performing a single Kempe switch: %d\n", maximum_different_pairs_reachable);
                print_colouring_from_vector(stderr, graph, factorisation_list[colouring_index]->factorisation);
                fprintf(stderr, "\n");
            }

            if(classes_per_graph) {
                fprintf(stderr, "%d equivalence class%s\n", uf_1factors->set_count,
                        uf_1factors->set_count == 1 ? "" : "es");

                int class_counter = 0;
                for (int i = 0; i < uf_1factors->size; ++i) {
                    if (unionfind_find_root(uf_1factors, i) == i) {
                        class_counter++;
                        fprintf(stderr, "Class %d:\n", class_counter);
                        fprintf(stderr, "  size: %d\n", uf_1factors->tree_size[i]);
                        int perfect_pair_count[k * (k - 1) / 2 + 1];
                        for (int j = 0; j < k * (k - 1) / 2 + 1; ++j) {
                            perfect_pair_count[j] = 0;
                        }
                        for (int j = 0; j < uf_1factors->size; ++j) {
                            if (unionfind_find_root(uf_1factors, j) == i)
                                perfect_pair_count[factorisation_list[j]->perfect_pair_count]++;
                        }
                        fprintf(stderr, "  perfect pairs: ");
                        for (int j = 0; j < k * (k - 1) / 2 + 1; ++j) {
                            if (perfect_pair_count[j]) {
                                fprintf(stderr, "%d(%d) ", j, perfect_pair_count[j]);
                            }
                        }
                        fprintf(stderr, "\n");
                    }
                }
                fprintf(stderr, "\n");
            }
            if(min_equivalence_class_count == -1 || uf_1factors->set_count < min_equivalence_class_count){
                min_equivalence_class_count = uf_1factors->set_count;
            }
            unionfind_free(uf_1factors);
            free_factorisation_list_element(factorisation_tree);
            free(factorisation_list);
            factorisation_tree = NULL;
        }

        free_plane_graph(graph);
    }

    fprintf(stderr, "Overview over all graphs\n========================\n\n");

    fprintf(stderr, "Read %d graph%s.\n\n", graph_count, graph_count==1 ? "" : "s");

    if(compute_equivalence_classes){
        fprintf(stderr, "Minimum number of Kempe equivalence classes over all graphs: %d\n\n", min_equivalence_class_count);

        fprintf(stderr, "Overview of changes in number of perfect pairs for single Kempe switches:\n");
        for(int i = 0; i < k*(k-1)/2+1; i++){
            fprintf(stderr, "%2d: ", i);
            for(int j = 0; j < k*(k-1)/2+1; j++){
                if(single_switch_attainable[i][j]){
                    fprintf(stderr, "%2d ", j);
                } else {
                    fprintf(stderr, "   ");
                }
            }
            fprintf(stderr, "\n");
        }

        fprintf(stderr, "\nMaximum different pairs reachable from one colouring with the specified number of perfect pairs:\n");
        for(int i = 0; i < k*(k-1)/2+1; i++){
            fprintf(stderr, "%2d: ", i);
            if(global_maximum_from_single_switch_per_count[i] < 0){
                fprintf(stderr, "_\n");
            } else {
                fprintf(stderr, "%d\n", global_maximum_from_single_switch_per_count[i]);
            }
        }
    }
}