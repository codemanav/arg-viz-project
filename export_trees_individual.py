import json, os

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
