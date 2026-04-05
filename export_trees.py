from ARGSimulator import sim_ARG
import json
import os

# Simulate a small ARG
ts = sim_ARG(sample_size=4, length=1000, plot=False)

def export_tskit_tables(ts, output_folder="output"):
    """
    Export nodes and edges from a tskit TreeSequence into JSON files
    suitable for visualization (D3.js / Cytoscape.js).
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
            "right": edge.right
        }
        for i, edge in enumerate(tables.edges)
    ]

    data = {"nodes": nodes, "edges": edges}

    with open(os.path.join(output_folder, "arg.json"), "w") as f:
        json.dump(data, f, indent=2)

    # Also export local trees as Newick files for phylotree.js
    for tree in ts.trees():
        start, end = map(int, tree.interval)
        with open(os.path.join(output_folder, f"tree_{start}_{end}.nwk"), "w") as f:
            f.write(tree.newick())

    print(f"✅ Exported ARG data to '{output_folder}/arg.json' and local trees as Newick files.")



def export_trees_manifest(ts, out_path="trees.json"):
    """
    Produces a single JSON file:
    [
      {"start": 0, "end": 198, "newick": "(...);"},
      {"start": 198, "end": 1000, "newick": "(...);"}
    ]
    """
    trees = []
    for tree in ts.trees():
        start, end = tree.interval
        trees.append({
            "start": float(start),
            "end": float(end),
            "newick": tree.newick()
        })
    with open(out_path, "w") as f:
        json.dump(trees, f, indent=2)
    print(f"Wrote {out_path} with {len(trees)} trees.")

# Export trees in Newick format
with open("local_trees.nwk", "w") as f:
    for tr_num, tree in enumerate(ts.trees()):
        newick_str = tree.newick()
        f.write(newick_str + ";\n")  # each tree separated by ;
        print(f"Tree {tr_num}: {newick_str}")

export_tskit_tables(ts)
export_trees_manifest(ts)