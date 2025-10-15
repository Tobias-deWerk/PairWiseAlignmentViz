# PairWiseAlignmentViz
A two-stream loop-and-bump visualization engine for the pairwise alignment of two sequences.

## Usage

```
python align_viz.py <input.fasta> <width> <height> <dpi> <output.svg>
```

The visualization depends on `matplotlib`. Install it in your active environment
before invoking the script (for example, `pip install matplotlib`).

Additional tuning parameters can be supplied through optional flags:

- `--min-gap-size`: gap length threshold used by weak-segment heuristics and gap summarization (default: 10).
- `--block-size`: maximum block length evaluated as a whole before switching to the sliding window heuristic (default: 100).
- `--min-sequence-identity`: minimum identity required to avoid flagging a segment as weakly aligned (default: 0.7).
- `--window-size`: window length for scanning long regions for weak alignment (default: 20).
- `--tick-interval`: spacing (in bp) for annotating local coordinate tick marks along each stream (default: 10,000; set to 0 to disable).
- `--backbone-gap`: vertical distance between the query and reference backbones (default: 1.0; lower values bring the streams closer together).
- `--backbone-thickness`: line width used for the backbone paths (default: 2.0).
- `--bump-scale`: multiplier applied to weak-alignment bump heights (default: 1.0).
- `--mismatch-line-width`: controls the thickness of mismatch ladder rungs (default: 1.2).
- `--gap-max-height`: maximum amplitude reached by gap beak glyphs (default: 0.8 data units).
- `--gap-width`: horizontal width assigned to each gap column on the global axis (default: 0.0 nucleotides).

All gaps are rendered as tapered beaks that originate at the backbone and arc
outward before rejoining the baseline. The glyph width and height scale with the
gap length up to a capped maximum (customizable through `--gap-max-height`) so
extremely long insertions remain legible without overwhelming nearby features.
If you would like the global axis to reflect gap length directly, set
`--gap-width` to a positive value to allocate horizontal space for each gap
column.

For troubleshooting, the tool also emits two TSV files alongside the requested
visualization output: `<output>_query_stream.tsv` and
`<output>_reference_stream.tsv`. Each row records the column index, global and
per-stream coordinates, local nucleotide positions, feature flags, and (when
applicable) the gap run length. This makes it straightforward to audit how the
gap glyphs were positioned if you need to refine their presentation.

The script accepts any pairwise alignment FASTA that contains exactly two sequences of equal length (including gap characters) such as MAFFT pairwise outputs.
