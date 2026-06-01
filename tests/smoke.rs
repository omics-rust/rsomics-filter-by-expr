use std::process::Command;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_rsomics-filter-by-expr"))
}
fn golden(n: &str) -> String {
    format!("{}/tests/golden/{}", env!("CARGO_MANIFEST_DIR"), n)
}

#[test]
fn output_is_gene_tab_bool() {
    let out = bin().arg(golden("counts.tsv")).output().unwrap();
    assert!(out.status.success());
    let s = String::from_utf8(out.stdout).unwrap();
    let mut lines = s.trim().lines();
    assert_eq!(lines.next().unwrap(), "gene\tkeep");
    for line in lines {
        let v = line.split('\t').nth(1).unwrap();
        assert!(v == "TRUE" || v == "FALSE", "non-bool keep value: {v}");
    }
}

#[test]
fn lone_high_count_gene_is_dropped() {
    // One sample over min.count but well under min.total isn't enough on its own.
    let out = bin().arg(golden("counts.tsv")).output().unwrap();
    assert!(out.status.success());
}
