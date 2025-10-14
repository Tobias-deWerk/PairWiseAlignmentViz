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
- `--backbone-gap`: vertical distance between the query and reference backbones (default: 1.0; lower values bring the streams closer together).
- `--backbone-thickness`: line width used for the backbone paths (default: 2.0).
- `--bump-scale`: multiplier applied to weak-alignment bump heights (default: 1.0).
- `--gap-circle-scale`: converts gap lengths (in bp) to the rendered circle diameter (default: 0.02; smaller values tighten the balloons).
- `--mismatch-line-width`: controls the thickness of mismatch ladder rungs (default: 1.2).

Large gaps are now rendered as balloon glyphs: each indel spawns a circle whose
diameter scales with the gap length and whose center sits at half the scaled
height above (query insertions) or below (reference insertions) the backbone.
The circles are jittered horizontally when necessary so they do not overlap, and
a curved stem ties each balloon back to the backbone at the precise gap
position. Every balloon receives a `+<length>` label so the inserted bases are
obvious at a glance. Shorter gaps (below `--min-gap-size`) continue to render as
compact beaks that hug the backbone.

For troubleshooting, the tool also emits two TSV files alongside the requested
visualization output: `<output>_query_stream.tsv` and
`<output>_reference_stream.tsv`. Each row records the column index, global and
per-stream coordinates, local nucleotide positions, feature flags, and (when
applicable) the gap circle geometry (anchor position, center, radius, and
applied jitter). This makes it straightforward to audit how the balloons were
scaled and displaced if you need to refine their presentation.

The script accepts any pairwise alignment FASTA that contains exactly two sequences of equal length (including gap characters) such as MAFFT pairwise outputs.
