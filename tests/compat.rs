use std::process::Command;

fn ours() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_BIN_EXE_rsomics-filter-by-expr"))
}

fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

fn rows(s: &str) -> Vec<(String, String)> {
    s.trim()
        .lines()
        .map(|l| {
            let mut f = l.split('\t');
            (
                f.next().unwrap().to_string(),
                f.next().unwrap_or("").to_string(),
            )
        })
        .collect()
}

// keep is boolean -> exact diff, no epsilon.
fn diff_exact(mine: &str, theirs: &str) {
    let a = rows(mine);
    let b = rows(theirs);
    assert_eq!(a.len(), b.len(), "row count mismatch");
    assert_eq!(a[0], b[0], "header mismatch");
    for (r, (x, y)) in a[1..].iter().zip(&b[1..]).enumerate() {
        assert_eq!(x.0, y.0, "row {r} gene id mismatch");
        assert_eq!(
            x.1, y.1,
            "row {r} gene {} keep mismatch: ours={} oracle={}",
            x.0, x.1, y.1
        );
    }
}

fn run_ours(extra: &[&str]) -> String {
    let out = Command::new(ours())
        .arg(golden("counts.tsv"))
        .args(extra)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn default_matches_golden() {
    diff_exact(
        &run_ours(&[]),
        &std::fs::read_to_string(golden("golden_keep.tsv")).unwrap(),
    );
}

#[test]
fn grouped_matches_golden() {
    diff_exact(
        &run_ours(&["--group", "a,a,a,b,b,b"]),
        &std::fs::read_to_string(golden("golden_keep_group.tsv")).unwrap(),
    );
}

// Live differential vs the edgeR R oracle. Loud-skips when conda/Rscript absent.
#[test]
fn matches_edger_oracle_live() {
    let Some(rscript) = conda_rscript() else {
        eprintln!("SKIP matches_edger_oracle_live: no `conda run -n r-bioc Rscript`");
        return;
    };
    let oracle = format!(
        "{}/tests/filter_by_expr_oracle.R",
        env!("CARGO_MANIFEST_DIR")
    );

    diff_exact(&run_ours(&[]), &run_oracle(&rscript, &oracle, &[]));
    diff_exact(
        &run_ours(&["--group", "a,a,a,b,b,b"]),
        &run_oracle(&rscript, &oracle, &["a,a,a,b,b,b"]),
    );
}

fn run_oracle(rscript: &[String], script: &str, extra: &[&str]) -> String {
    let mut cmd = Command::new(&rscript[0]);
    cmd.args(&rscript[1..])
        .arg(script)
        .arg(golden("counts.tsv"));
    for e in extra {
        cmd.arg(e);
    }
    let out = cmd.output().unwrap();
    assert!(
        out.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn conda_rscript() -> Option<Vec<String>> {
    let home = std::env::var("HOME").unwrap_or_default();
    let candidates: Vec<Vec<String>> = vec![
        vec![
            "conda".into(),
            "run".into(),
            "-n".into(),
            "r-bioc".into(),
            "Rscript".into(),
        ],
        vec![
            format!("{home}/miniconda3/bin/conda"),
            "run".into(),
            "-n".into(),
            "r-bioc".into(),
            "Rscript".into(),
        ],
    ];
    for c in candidates {
        let ok = Command::new(&c[0])
            .args(&c[1..])
            .arg("-e")
            .arg("cat('ok')")
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(c);
        }
    }
    None
}
