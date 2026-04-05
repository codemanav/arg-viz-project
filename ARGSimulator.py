from Espalier import MAF
import msprime
import dendropy
import json
import os

def sim_ARG(sample_size=10, Ne=100, length=1e3, recombination_rate=5e-6,
            min_breakpoints=2, max_breakpoints=3, plot=False):
    """
    Simulate ARGs with local trees in a tree sequence using msprime.

    Parameters:
        sample_size (int): Number of tips/samples in ARG
        Ne (float): Haploid effective population size
        length (int): Genome length
        recombination_rate (float): recombination rate per site per lineage
        min_breakpoints (int): minimum number of allowed breakpoints
        max_breakpoints (int): maximum number of allowed breakpoints
        plot (boolean): display ARG as local trees and ts tables

    Returns:
        ts (tskit.TreeSequence): TreeSequence representation of simulated ARG
    """
    breaks = 0
    while breaks < min_breakpoints or breaks > max_breakpoints:
        ts = msprime.simulate(
            sample_size=sample_size, Ne=Ne, length=length,
            recombination_rate=recombination_rate, record_full_arg=True
        )
        breaks = len(ts.breakpoints(as_array=True)) - 2

    for tree in ts.trees():
        if plot:
            print("-" * 20)
            print("tree {}: interval = {}".format(tree.index, tree.interval))
            print(tree.draw(format="unicode"))

    topo_changes = count_topo_changes(ts)

    if plot:
        print(ts.tables.nodes)
        print()
        print(ts.tables.edges)
        print()
        print("Recombination breakpoints: " + str(breaks))
        print()
        print("Topological changes in tree sequence: " + str(topo_changes))

    export_tskit_tables(ts)
    return ts


def export_tskit_tables(ts, output_folder="output"):
    """
    Export nodes, edges, and local trees from a tskit TreeSequence
    into JSON files consumed by the D3.js visualizations.

    Produces:
        arg.json   — full ARG node/edge tables
        trees.json — marginal local trees with Newick strings
    """
    os.makedirs(output_folder, exist_ok=True)
    tables = ts.dump_tables()

    nodes = [
        {"id": i, "time": node.time, "flags": int(node.flags)}
        for i, node in enumerate(tables.nodes)
    ]

    edges = [
        {
            "id": i,
            "parent": int(edge.parent),
            "child": int(edge.child),
            "left": edge.left,
            "right": edge.right,
        }
        for i, edge in enumerate(tables.edges)
    ]

    with open(os.path.join(output_folder, "arg.json"), "w") as f:
        json.dump({"nodes": nodes, "edges": edges}, f, indent=2)

    trees_data = []
    for tree in ts.trees():
        start, end = tree.interval
        trees_data.append({
            "start": float(start), "end": float(end),
            "newick": tree.newick(),
        })

    with open(os.path.join(output_folder, "trees.json"), "w") as f:
        json.dump(trees_data, f, indent=2)

    print(f"✅ Exported to '{output_folder}/': arg.json, trees.json")


def count_topo_changes(ts):
    """
    Counts topological changes between adjacent local trees
    using Espalier's MAF discordance test.

    Parameters:
        ts (tskit.TreeSequence): Tree sequence

    Returns:
        changes (int): number of topological changes
    """
    changes = 0
    prev_tree = None
    for tr_num, tree in enumerate(ts.trees()):
        if tr_num > 0:
            taxa = dendropy.TaxonNamespace()
            ref = dendropy.Tree.get(
                data=prev_tree, schema="newick",
                rooting="default-rooted", taxon_namespace=taxa
            )
            alt = dendropy.Tree.get(
                data=tree.newick(), schema="newick",
                rooting="default-rooted", taxon_namespace=taxa
            )
            if MAF.test_discordance(ref, alt):
                changes += 1
        prev_tree = tree.newick()

    return changes
