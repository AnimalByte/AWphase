use std::fs;

use anyhow::{Context, Result};

use crate::types::{ReadObservation, WindowId};

pub fn load_reads_json(path: &str, window: &WindowId) -> Result<Vec<ReadObservation>> {
    let text =
        fs::read_to_string(path).with_context(|| format!("failed to read reads file: {path}"))?;
    let all: Vec<ReadObservation> = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse reads JSON: {path}"))?;

    Ok(all
        .into_iter()
        .filter(|r| r.site_pos >= window.start && r.site_pos <= window.end)
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str) -> PathBuf {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("awphase_{name}_{ts}.json"))
    }

    #[test]
    fn loads_and_filters_reads_by_window() {
        let path = temp_path("reads");
        let json = r#"
[
  {"read_id":"r1","site_pos":101,"allele":1,"baseq":30,"mapq":60},
  {"read_id":"r2","site_pos":251,"allele":-1,"baseq":30,"mapq":60},
  {"read_id":"r3","site_pos":9999,"allele":1,"baseq":30,"mapq":60}
]
"#;
        std::fs::write(&path, json).unwrap();

        let window = WindowId {
            chrom: "chr20".to_string(),
            start: 100,
            end: 300,
        };

        let reads = load_reads_json(path.to_str().unwrap(), &window).unwrap();
        assert_eq!(reads.len(), 2);
        assert_eq!(reads[0].read_id, "r1");
        assert_eq!(reads[1].read_id, "r2");

        let _ = std::fs::remove_file(path);
    }
}
