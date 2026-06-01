# rsomics-filter-by-expr

Boolean per-gene expression filter for an RNA-seq count matrix — a Rust
reimplementation of edgeR's `filterByExpr`.

Reads a gene x sample count matrix (TSV: first column gene id, header row of
sample names) and writes `gene<TAB>keep` where `keep` is `TRUE`/`FALSE`.

```
rsomics-filter-by-expr counts.tsv -o keep.tsv
rsomics-filter-by-expr counts.tsv --group ctrl,ctrl,trt,trt -o keep.tsv
```

A gene is kept when it clears both:

- a CPM cutoff (`min.count` scaled to a CPM via the median library size) in at
  least the required number of samples, and
- a minimum summed count across all samples (`min.total.count`).

The required sample count is the smallest group size (all samples when no group
is given), tapered past `large.n` by `min.prop`. Library sizes default to column
sums; pass `--lib-size` to supply your own (e.g. raw size times a TMM factor).

Defaults match edgeR: `min.count=10`, `min.total.count=15`, `large.n=10`,
`min.prop=0.7`.

## Origin

Independent Rust reimplementation of edgeR's `filterByExpr.default`
(Bioconductor). edgeR's source is GPL-licensed; this crate reads only the
public function definition and reproduces the documented default keep rule:

- CPM cutoff `= min.count / median(lib.size) * 1e6`
- `keep.CPM = rowSums(CPM >= cutoff) >= (MinSampleSize - tol)`
- `keep.TotalCount = rowSums(y) >= min.total.count - tol`
- `MinSampleSize` = smallest non-empty group, capped at
  `large.n + (n - large.n) * min.prop`

Correctness is validated value-exact against the edgeR R binary
(`tests/compat.rs` + `tests/filter_by_expr_oracle.R`); committed goldens let CI
validate without R installed.

Reference: Chen Y, Lun ATL, Smyth GK. *F1000Research* 2016 (edgeR);
McCarthy DJ, Chen Y, Smyth GK. *Nucleic Acids Research* 2012.

License: MIT OR Apache-2.0.
Upstream credit: edgeR https://bioconductor.org/packages/edgeR (GPL >= 2).
