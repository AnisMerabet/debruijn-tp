#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
matplotlib.use("Agg")

__author__ = "Merabet"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Merabet"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Merabet"
__email__ = "merabet9anis@gmail.com"
__status__ = "Developpement"

def isfile(path): # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file
    
    :raises ArgumentTypeError: If file doesn't exist
    
    :return: (str) Path 
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="k-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file (default contigs.fasta)")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as an image (png)")
    return parser.parse_args()


def read_fastq(fastq_file):
    """Extract reads from fastq files.

    :param fastq_file: (str) Path to the fastq file.
    :return: A generator object that iterate the read sequences. 
    """
    with open(fastq_file, 'r') as fastq_file:
        while True:
            try:
                # Skip header line
                header = next(fastq_file)
                # Read sequence line
                sequence = next(fastq_file).strip()
                # Skip the quality line
                quality = next(fastq_file)
                # Skip the plus line
                plus = next(fastq_file)
                yield sequence
            except StopIteration:
                break


def cut_kmer(read, kmer_size):
    """Cut read into kmers of size kmer_size.
    
    :param read: (str) Sequence of a read.
    :return: A generator object that iterate the kmers of of size kmer_size.
    """
    for i in range(0, len(read) - kmer_size + 1):
        yield read[i:i + kmer_size]


def build_kmer_dict(fastq_file, kmer_size):
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}  # Initialize an empty dictionary to store k-mer occurrences
    
    # Iterate through each read sequence in the fastq file
    for read_sequence in read_fastq(fastq_file):
        # Generate k-mers from the read sequence
        for kmer in cut_kmer(read_sequence, kmer_size):
            kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    
    return kmer_dict

def build_graph(kmer_dict):
    """Build the debruijn graph

    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    # Create a directed graph
    debruijn_graph = nx.DiGraph()
    
    kmer_items = list(kmer_dict.items())
    
    for i, (kmer, count) in enumerate(kmer_items):
        debruijn_graph.add_edge(kmer[:-1], kmer[1:], weight=count)
    
    return debruijn_graph


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        for ntd in range(len(path) - 1):

            graph.remove_edge(path[ntd], path[ntd+1])

    if delete_entry_node and delete_sink_node:
        # Remove all nodes in the path
        graph.remove_nodes_from([i[0] for i in path_list])
        graph.remove_nodes_from([i[-1] for i in path_list])

    graph.remove_nodes_from(list(nx.isolates(graph)))

    return graph

def select_best_path(graph, path_list, path_length, weight_avg_list, 
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path 
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    if len(path_list) == 1:
        # If there's only one path, return the original graph
        return graph.copy()
    
    # Calculate standard deviation of average weights and path lengths
    weight_stddev = statistics.stdev(weight_avg_list)
    length_stddev = statistics.stdev(path_length)
    
    if weight_stddev > 0:
        # If weights vary, select the path with the highest weight
        best_path_index = weight_avg_list.index(max(weight_avg_list))
    elif length_stddev > 0:
        # If weights are equal but lengths vary, select the longest path
        best_path_index = path_length.index(max(path_length))
    else:
        # If weights and lengths are equal, choose randomly
        best_path_index = random.randint(0, len(path_list) - 1)
    
    # Remove other paths from the graph based on the best path index
    for i, path in enumerate(path_list):
        if i != best_path_index:
            modified_graph = remove_paths(graph, [path], delete_entry_node, delete_sink_node)
    
    return modified_graph


def path_average_weight(graph, path):
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph 
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    all_paths = list(nx.all_simple_paths(graph, ancestor_node, descendant_node))
    all_weigths = [path_average_weight(graph, path) for path in all_paths]
    all_length = [len(path)-1 for path in all_paths]
    graph = select_best_path(graph, all_paths, all_length, all_weigths)
    return graph

def simplify_bubbles(graph):
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    bubble_detected = False
    bubble_ancestor = None
    bubble_descendant = None

    for node in graph:
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for i, node_i in enumerate(predecessors):
                for node_j in predecessors[i + 1:]:
                    common_ancestor = nx.lowest_common_ancestor(graph, node_i, node_j)
                    if common_ancestor is not None:
                        bubble_ancestor = common_ancestor
                        bubble_descendant = node
                        bubble_detected = True
                        break
                if bubble_detected:
                    break

    if bubble_detected:
        graph = simplify_bubbles(solve_bubble(graph, bubble_ancestor, bubble_descendant))

    return graph

def solve_entry_tips(graph, starting_nodes):
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    while True:
        entry_tips = [node for node in graph if len(list(graph.predecessors(node))) > 1]

        if not entry_tips:
            break

        for node in entry_tips:
            paths = [list(nx.all_simple_paths(graph, node_start_i, node)) for node_start_i in starting_nodes]
            paths = [path[0] for path in paths if len(path) > 0]
            lengths = [len(path) - 1 for path in paths]
            weights = [path_average_weight(graph, path) if lengths[i] > 1 else graph.get_edge_data(*path)["weight"]
                       for i, path in enumerate(paths)]

            graph = select_best_path(graph, paths, lengths, weights,
                                     delete_entry_node=True, delete_sink_node=False)

    return graph

def solve_out_tips(graph, ending_nodes):
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    while True:
        found_tip = False
        for node in graph:
            node_success = list(graph.successors(node))
            if len(node_success) > 1:
                paths = [list(nx.all_simple_paths(graph, node, node_end_i))\
                         for node_end_i in ending_nodes]
                paths = [path[0] for path in paths if len(path) > 0]
                lengths = [len(path) - 1 for path in paths]
                weights = [path_average_weight(graph, path) if lengths[i] > 1 else \
                           graph.get_edge_data(*path)["weight"]
                           for i, path in enumerate(paths)]

                graph = select_best_path(graph, paths, lengths, weights, 
                                         delete_entry_node=False, delete_sink_node=True)
                found_tip = True
                break

        if not found_tip:
            break

    return graph

def get_starting_nodes(graph):
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    starting_nodes = [node for node in graph.nodes() if len(list(graph.predecessors(node))) == 0]
    return starting_nodes

def get_sink_nodes(graph):
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    sink_nodes = [node for node in graph.nodes() if len(list(graph.successors(node))) == 0]
    return sink_nodes

def get_contigs(graph, starting_nodes, ending_nodes):
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object 
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []
    
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            # Check if there is a path from a starting node to an ending node
            if nx.has_path(graph, start_node, end_node):
                # Find all simple paths (contigs) between the starting and ending nodes
                for path in nx.all_simple_paths(graph, start_node, end_node):
                    contig = path[0]
                    for i in range(1, len(path)):
                        contig += path[i][-1]  # Add the last character of each node to the contig
                    contigs.append((contig, len(contig)))

    return contigs

def save_contigs(contigs_list, output_file):
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (str) Path to the output file
    """
    with open(output_file, 'w') as file:
        for i, (contig, length) in enumerate(contigs_list):
            header = f">contig_{i} len={length}\n"
            formatted_sequence = textwrap.fill(contig, width=80)
            file.write(header + formatted_sequence + '\n')


def draw_graph(graph, graphimg_file): # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (str) Path to the output file
    """                                   
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]
    #print(elarge)
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    #print(elarge)
    # Draw the graph with networkx
    #pos=nx.spring_layout(graph)
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5, 
                           edge_color='b', style='dashed')
    #nx.draw_networkx(graph, pos, node_size=10, with_labels=False)
    # save image
    plt.savefig(graphimg_file)


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    # graphe
    args = get_arguments()
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    graph = simplify_bubbles(graph)
    graph = solve_entry_tips(graph, get_starting_nodes(graph))
    graph = solve_out_tips(graph, get_sink_nodes(graph))
    contigs = get_contigs(graph, get_starting_nodes(graph), get_sink_nodes(graph))
    save_contigs(contigs, args.output_file)

    # Plot the graph
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)

if __name__ == '__main__': # pragma: no cover
    main()
