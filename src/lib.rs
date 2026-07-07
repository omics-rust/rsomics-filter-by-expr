use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub struct Matrix {
    pub gene_col: String,
    pub genes: Vec<String>,
    pub counts: Vec<f64>,
    pub n_samples: usize,
}

impl Matrix {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();

        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty count matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let mut hdr = header.split('\t');
        let gene_col = hdr
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("header without a gene column".into()))?
            .to_string();
        let n_samples = hdr.count();
        if n_samples == 0 {
            return Err(RsomicsError::InvalidInput(
                "header has no sample columns".into(),
            ));
        }

        let mut genes = Vec::new();
        let mut counts = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let gene = fields
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("row without a gene id".into()))?;
            genes.push(gene.to_string());
            let before = counts.len();
            for f in fields {
                let v = f.parse::<f64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric count '{f}' for gene {gene}"))
                })?;
                if !v.is_finite() {
                    return Err(RsomicsError::InvalidInput(format!(
                        "non-finite count '{f}' for gene {gene}"
                    )));
                }
                if v < 0.0 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "gene {gene}: negative counts not allowed"
                    )));
                }
                counts.push(v);
            }
            if counts.len() - before != n_samples {
                return Err(RsomicsError::InvalidInput(format!(
                    "gene {gene}: {} values, header has {n_samples} samples",
                    counts.len() - before
                )));
            }
        }
        Ok(Self {
            gene_col,
            genes,
            counts,
            n_samples,
        })
    }
}

fn load_lib_sizes(path: &Path, n_samples: usize) -> Result<Vec<f64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut sizes = Vec::with_capacity(n_samples);
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let val = line.rsplit('\t').next().unwrap_or(line);
        sizes
            .push(val.parse::<f64>().map_err(|_| {
                RsomicsError::InvalidInput(format!("non-numeric lib size '{val}'"))
            })?);
    }
    if sizes.len() != n_samples {
        return Err(RsomicsError::InvalidInput(format!(
            "{} lib sizes for {n_samples} samples",
            sizes.len()
        )));
    }
    Ok(sizes)
}

fn column_sums(m: &Matrix) -> Vec<f64> {
    let mut sizes = vec![0.0f64; m.n_samples];
    for row in m.counts.chunks_exact(m.n_samples) {
        for (s, &c) in sizes.iter_mut().zip(row) {
            *s += c;
        }
    }
    sizes
}

// R's stats::median: sort, then mean of the two central order statistics for even n.
fn median(values: &[f64]) -> f64 {
    let mut v = values.to_vec();
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

pub struct Defaults {
    pub min_count: f64,
    pub min_total_count: f64,
    pub large_n: f64,
    pub min_prop: f64,
}

impl Default for Defaults {
    fn default() -> Self {
        Self {
            min_count: 10.0,
            min_total_count: 15.0,
            large_n: 10.0,
            min_prop: 0.7,
        }
    }
}

pub struct FilterOpts {
    pub defaults: Defaults,
    pub group: Option<Vec<usize>>,
}

// edgeR filterByExpr.default with no design: the smallest group decides how many
// samples must clear the CPM cutoff, capped at large.n via the min.prop taper.
fn min_sample_size(opts: &FilterOpts, n_samples: usize) -> f64 {
    let d = &opts.defaults;
    let mut mss = match &opts.group {
        Some(g) => {
            let n_groups = g.iter().copied().max().map(|m| m + 1).unwrap_or(0);
            let mut counts = vec![0usize; n_groups];
            for &gi in g {
                counts[gi] += 1;
            }
            counts.into_iter().filter(|&c| c > 0).min().unwrap_or(0) as f64
        }
        None => n_samples as f64,
    };
    if mss > d.large_n {
        mss = d.large_n + (mss - d.large_n) * d.min_prop;
    }
    mss
}

pub fn filter_by_expr(
    counts_path: &Path,
    lib_size_path: Option<&Path>,
    opts: &FilterOpts,
    output: &mut dyn Write,
) -> Result<(u64, u64)> {
    let m = Matrix::load(counts_path)?;
    let lib_size = match lib_size_path {
        Some(p) => load_lib_sizes(p, m.n_samples)?,
        None => column_sums(&m),
    };
    if lib_size.iter().any(|&l| l <= 0.0) {
        return Err(RsomicsError::InvalidInput(
            "library sizes should be greater than zero".into(),
        ));
    }

    let d = &opts.defaults;
    let mss = min_sample_size(opts, m.n_samples);
    let median_lib = median(&lib_size);
    let cpm_cutoff = d.min_count / median_lib * 1e6;
    // edgeR cpm.default: count / (lib_size / 1e6), not count * (1e6 / lib_size);
    // the two differ by 1 ULP and flip keep at a boundary CPM == cutoff.
    let cpm_denom: Vec<f64> = lib_size.iter().map(|&l| l / 1e6).collect();

    const TOL: f64 = 1e-14;
    let n_required = mss - TOL;
    let total_cutoff = d.min_total_count - TOL;

    let mut out = BufWriter::new(output);
    out.write_all(m.gene_col.as_bytes())
        .map_err(RsomicsError::Io)?;
    out.write_all(b"\tkeep\n").map_err(RsomicsError::Io)?;

    let mut kept = 0u64;
    for (gene, row) in m.genes.iter().zip(m.counts.chunks_exact(m.n_samples)) {
        let mut above = 0usize;
        let mut total = 0.0f64;
        for (&c, &denom) in row.iter().zip(&cpm_denom) {
            total += c;
            if c / denom >= cpm_cutoff {
                above += 1;
            }
        }
        let keep = (above as f64) >= n_required && total >= total_cutoff;
        if keep {
            kept += 1;
        }
        out.write_all(gene.as_bytes()).map_err(RsomicsError::Io)?;
        out.write_all(if keep { b"\tTRUE\n" } else { b"\tFALSE\n" })
            .map_err(RsomicsError::Io)?;
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok((m.genes.len() as u64, kept))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_even_and_odd() {
        assert_eq!(median(&[3.0, 1.0, 2.0]), 2.0);
        assert_eq!(median(&[4.0, 1.0, 3.0, 2.0]), 2.5);
    }

    #[test]
    fn min_sample_size_one_group_is_all_samples_tapered() {
        let opts = FilterOpts {
            defaults: Defaults::default(),
            group: None,
        };
        // 4 samples <= large.n=10, no taper.
        assert_eq!(min_sample_size(&opts, 4), 4.0);
        // 20 samples > 10: 10 + (20-10)*0.7 = 17.
        assert_eq!(min_sample_size(&opts, 20), 17.0);
    }

    #[test]
    fn smallest_group_decides() {
        let opts = FilterOpts {
            defaults: Defaults::default(),
            group: Some(vec![0, 0, 0, 1, 1]),
        };
        assert_eq!(min_sample_size(&opts, 5), 2.0);
    }
}
