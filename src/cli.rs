use std::collections::BTreeMap;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Section};

use rsomics_filter_by_expr::{Defaults, FilterOpts, filter_by_expr};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-filter-by-expr", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[arg(long)]
    group: Option<String>,
    #[arg(long)]
    lib_size: Option<PathBuf>,
    #[arg(long, default_value_t = 10.0)]
    min_count: f64,
    #[arg(long, default_value_t = 15.0)]
    min_total_count: f64,
    #[arg(long, default_value_t = 10.0)]
    large_n: f64,
    #[arg(long, default_value_t = 0.7)]
    min_prop: f64,
    #[command(flatten)]
    pub common: CommonFlags,
}

fn parse_group(spec: &str) -> Vec<usize> {
    let mut ids = BTreeMap::new();
    let mut next = 0usize;
    spec.split(',')
        .map(|label| {
            *ids.entry(label.to_string()).or_insert_with(|| {
                let id = next;
                next += 1;
                id
            })
        })
        .collect()
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let mut out: Box<dyn std::io::Write> = if self.output == "-" && self.common.json {
            Box::new(std::io::sink())
        } else if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        let opts = FilterOpts {
            defaults: Defaults {
                min_count: self.min_count,
                min_total_count: self.min_total_count,
                large_n: self.large_n,
                min_prop: self.min_prop,
            },
            group: self.group.as_deref().map(parse_group),
        };
        let (total, kept) =
            filter_by_expr(&self.counts, self.lib_size.as_deref(), &opts, &mut out)?;
        if !self.common.quiet {
            eprintln!("{kept}/{total} genes kept");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Boolean per-gene expression filter for a count matrix (edgeR filterByExpr).",
    origin: None,
    usage_lines: &["<counts.tsv> [--group a,a,b,b] [--lib-size f.tsv] [-o keep.tsv]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "group",
                aliases: &[],
                value: Some("<labels>"),
                type_hint: Some("String"),
                required: false,
                default: None,
                description: "Comma-separated group label per sample; the smallest group sets the sample threshold.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "lib-size",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-sample library sizes (e.g. raw size times TMM factor); defaults to column sums.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "min-count",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("10"),
                description: "Per-sample count behind the CPM cutoff.",
                why_default: Some("edgeR default min.count."),
            },
            FlagSpec {
                short: None,
                long: "min-total-count",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("15"),
                description: "Minimum summed count across all samples.",
                why_default: Some("edgeR default min.total.count."),
            },
            FlagSpec {
                short: None,
                long: "large-n",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("10"),
                description: "Group size above which the sample threshold tapers.",
                why_default: Some("edgeR default large.n."),
            },
            FlagSpec {
                short: None,
                long: "min-prop",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.7"),
                description: "Taper proportion applied past large.n.",
                why_default: Some("edgeR default min.prop."),
            },
        ],
    }],
    examples: &[
        Example {
            description: "Default filter, one group",
            command: "rsomics-filter-by-expr counts.tsv -o keep.tsv",
        },
        Example {
            description: "Two-group design",
            command: "rsomics-filter-by-expr counts.tsv --group ctrl,ctrl,trt,trt -o keep.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }

    #[test]
    fn group_labels_map_to_first_seen_order() {
        assert_eq!(parse_group("a,a,b,b,c"), vec![0, 0, 1, 1, 2]);
        assert_eq!(parse_group("trt,ctrl,trt"), vec![0, 1, 0]);
    }
}
