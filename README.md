# ARG Visualization

Web-based visualization of simulated **Ancestral Recombination Graphs (ARGs)** and their marginal local trees, targeting the figure style from [Rasmussen et al., PLOS Genetics 2014](https://doi.org/10.1371/journal.pgen.1004342).

## Project structure

```
ARG/
├── ARGSimulator.py              # Data generation (msprime → JSON)
├── output/
│   ├── arg.json                 # Full ARG: nodes + interval-annotated edges
│   ├── trees.json               # Marginal local trees (Newick per interval)
│   ├── arg_rasmussen.html       # Visualization 1 — Rasmussen-style figure
│   └── show_latest.html         # Visualization 2 — interactive dashboard
├── images/
│   ├── rasmussen style visual.png   # Screenshot of arg_rasmussen.html
│   └── Interactive Dasboard.png     # Screenshot of show_latest.html
└── README.md
```

## Requirements

Python 3.9+ with the following packages:

```
msprime
dendropy
Espalier
```

Install all at once:

```bash
pip install msprime dendropy Espalier
```

## Usage

### 1. Generate data

```python
from ARGSimulator import sim_ARG

ts = sim_ARG(sample_size=4, length=1000)
```

This writes `output/arg.json` and `output/trees.json` from the same simulation run.

`sim_ARG` parameters:

| Parameter            | Default  | Description                              |
|----------------------|----------|------------------------------------------|
| `sample_size`        | 10       | Number of sampled haplotypes             |
| `Ne`                 | 100      | Haploid effective population size        |
| `length`             | 1000     | Genome length in base pairs              |
| `recombination_rate` | 5e-6     | Per-site per-lineage recombination rate  |
| `min_breakpoints`    | 2        | Minimum recombination breakpoints        |
| `max_breakpoints`    | 3        | Maximum recombination breakpoints        |
| `plot`               | False    | Print ASCII trees and tskit tables       |

### 2. Serve the visualizations

Both HTML files load JSON via `fetch()`, so they must be served over HTTP:

```bash
cd output
python3 -m http.server 8765
```

Then open in a browser:

- **Rasmussen-style figure** — [http://localhost:8765/arg_rasmussen.html](http://localhost:8765/arg_rasmussen.html)
- **Interactive dashboard** — [http://localhost:8765/show_latest.html](http://localhost:8765/show_latest.html)

## Visualizations

### `arg_rasmussen.html` — Rasmussen-style figure (D3.js + SVG)

![Rasmussen-style ARG visualization](images/rasmussen%20style%20visual.png)

Replicates the two-panel layout from Rasmussen et al. 2014, Figure 1.

**Panel A — Full ARG.** Samples at the bottom, root at the top, time axis on the left (sqrt-scaled). Each genomic segment between recombination breakpoints is assigned a distinct color. Edges are drawn as orthogonal connectors (vertical then horizontal) with parallel colored lines showing which genomic intervals each edge spans. Genome color bars under each sample and a central genome ruler label the breakpoints.

**Panel B — Local trees.** One rooted tree per genomic segment, drawn side by side in the same color, with recombination labels (R) and arrows between them. A shared genome ruler at the bottom maps each tree to its interval.

Data consumed: **`arg.json`** only (breakpoints, segments, and local trees are derived from the edge intervals).

### `show_latest.html` — Interactive dashboard (Cytoscape + D3.js)

![Interactive dashboard](images/Interactive%20Dasboard.png)

Split-pane layout. Left panel renders the full ARG as a Cytoscape.js graph with dendrogram-style preset positions. Right panel draws the currently selected local tree as a D3 phylogram. A slider selects the genomic tree index; active edges are highlighted in blue on the ARG.

Data consumed: **`arg.json`** and **`trees.json`**.

## Data formats

### `arg.json`

```json
{
  "nodes": [
    { "id": 0, "time": 0.0, "flags": 1 }
  ],
  "edges": [
    { "id": 0, "parent": 4, "child": 0, "left": 0.0, "right": 1000.0 }
  ]
}
```

- **`flags`**: `1` = sample, `0` = coalescent, `131072` = recombination, `262144` = common ancestor event
- **`left` / `right`**: genomic interval over which the edge is active

### `trees.json`

```json
[
  { "start": 0.0, "end": 496.0, "newick": "((1:38,3:38):64,(0:75):27);" }
]
```

One entry per marginal tree, covering `[start, end)` of the genome.

## Pipeline overview

```
msprime.simulate(record_full_arg=True)
        │
        ▼
  tskit.TreeSequence
        │
        ├──► arg.json   (full DAG: nodes + interval-annotated edges)
        │
        └──► trees.json (marginal local trees as Newick strings)
                │
                ▼
   ┌────────────────────────────┐
   │  arg_rasmussen.html        │  ← reads arg.json
   │  show_latest.html          │  ← reads arg.json + trees.json
   └────────────────────────────┘
```
