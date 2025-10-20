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
- `--gap-max-height`: maximum amplitude reached by the gap glyphs (default: 0.8 data units).
- `--gap-width`: horizontal width assigned to every gap glyph (default: 0.0 nucleotides).
- `--gap-height-scale`: multiplier applied to gap amplitudes before capping at the maximum height (default: 0.04).
- `--indel-height-scale`: multiplier used when computing the height of short bird-beak indels (default: 0.04).
- `--gap-label-size`: font size for the `+x bp` annotations placed on gap glyphs (default: 8.0; pass `NA` or `NULL` to hide labels).
- `--query-annotation`: path to a per-stream gene annotation file for the query backbone. Each annotation file contains tab-delimited gene and feature records (see below).
- `--reference-annotation`: path to a per-stream gene annotation file for the reference backbone.
- `--annotation-label-size`: font size for the gene labels drawn alongside annotations (default: 10.0; pass `NA` or `NULL` to hide labels).
- `--annotation-thickness`: line width used for coding segments in gene annotations (default: 4.0).
- `--annotation-alpha`: alpha transparency applied to annotation overlays (default: 0.85).
- `--reference-annotation-color`: color used for reference annotations (default: `#984ea3`).
- `--query-annotation-color`: color used for query annotations (default: `#ff7f00`).
- `--annotation-label-jitter`: horizontal jitter amplitude for annotation labels to reduce overlap (default: 0.35; set to 0 to disable).
- `--annotation-max-layers`: maximum number of stacked annotation tracks drawn per stream (default: 3).
- `--annotation-spacing`: vertical distance (in data units) between stacked annotation tracks (default: 0.8; set to 0 to collapse them onto the backbone).

Insertions that meet or exceed `--min-gap-size` are rendered as smooth, fixed-width
Bezier arcs that peel away from the backbone and rejoin it with a mirrored slope.
Shorter indels continue to use zero-width bird-beak spikes so they remain crisp and
compact. The amplitude of long gaps is derived from `--gap-height-scale × length`,
while short indels use `--indel-height-scale × length`; both are clamped to
`--gap-max-height`. `--gap-width` controls the horizontal span
allocated to each large gap, while optional `+x bp` labels (suppressed when the
label size is `NA`) are jittered slightly along the x-axis to reduce overlap.
Set the width to zero to keep every gap anchored to a single global coordinate,
or increase it to open space for the Bezier curves and their annotations.

For troubleshooting, the tool also emits two TSV files alongside the requested
visualization output: `<output>_query_stream.tsv` and
`<output>_reference_stream.tsv`. Each row records the column index, global and
per-stream coordinates, local nucleotide positions, feature flags, and (when
applicable) the gap run length. This makes it straightforward to audit how the
gap glyphs were positioned if you need to refine their presentation.

The script accepts any pairwise alignment FASTA that contains exactly two sequences of equal length (including gap characters) such as MAFFT pairwise outputs.

## Gene annotations

When `--query-annotation` and/or `--reference-annotation` are supplied, the visualization overlays per-gene features on the respective stream. Annotation files follow a block-based format:

```
<gene_id>\t<gene_start>\t<gene_end>\t<direction>
<feature_name>\t<feature_start>\t<feature_end>
<feature_name>\t<feature_start>\t<feature_end>

<next_gene_id> ...
```

Gene coordinates are expressed in the local coordinate system for the stream (1-based, ignoring gaps). The direction column must be either `<` or `>` and is included in the rendered labels (e.g., `< gene1` or `gene2 >`). Feature rows describe named segments—commonly `UTR` and `exon` entries. Any gaps between features are automatically rendered as introns. Exons are drawn with a thick overlay, whereas UTRs and introns inherit the backbone width. Overlapping genes are automatically assigned to stacked annotation tracks (up to `--annotation-max-layers` per stream), each offset from the backbone by `--annotation-spacing`. Additional overlaps beyond the configured layer cap are omitted. Labels can be disabled by passing `--annotation-label-size NULL`, and their horizontal jitter can be tuned with `--annotation-label-jitter`.
