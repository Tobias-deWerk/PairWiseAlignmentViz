# PairWiseAlignmentViz
A two-stream loop-and-bump visualization engine for the pairwise alignment of two sequences.

## Usage

```
python align_viz.py <input.fasta> <width> <height> <dpi> <output.svg>
```

The visualization depends on `matplotlib`. Install it in your active environment
before invoking the script (for example, `pip install matplotlib`).

Additional tuning parameters can be supplied through optional flags:

- `--min-gap-size`: gap length threshold separating loop glyphs from bird beaks (default: 10).
- `--block-size`: maximum block length evaluated as a whole before switching to the sliding window heuristic (default: 100).
- `--min-sequence-identity`: minimum identity required to avoid flagging a segment as weakly aligned (default: 0.7).
- `--window-size`: window length for scanning long regions for weak alignment (default: 20).
- `--tick-interval`: spacing (in bp) for annotating local coordinate tick marks along each stream (default: 10,000; set to 0 to disable).

Gap loops now bloom to alternating sides so the strands do not collide with their
neighbors. Their shapes scale smoothly with the indel length while remaining
bounded, preventing the "pulled string" artefacts that previously dragged long
insertions across the figure. A `+<length>` annotation is placed at the apex of
each loop (for gaps that exceed `--min-gap-size`) so the amount of inserted
sequence remains immediately visible even though the rendered curve is compact.
The debug output (described below) also records the actual loop geometry metrics
for downstream inspection.

For troubleshooting, the tool also emits two TSV files alongside the requested
visualization output: `<output>_query_stream.tsv` and
`<output>_reference_stream.tsv`. Each row records the column index, global and
per-stream coordinates, local nucleotide positions, feature flags, and (when
applicable) the loop size, amplitude, width, displayed arc length, and blooming
direction. This makes it straightforward to audit how loop glyphs were scaled
and oriented if you need to refine their presentation.

The script accepts any pairwise alignment FASTA that contains exactly two sequences of equal length (including gap characters) such as MAFFT pairwise outputs.
