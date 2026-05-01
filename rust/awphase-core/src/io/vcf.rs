use std::fs;

use anyhow::{Context, Result};

use crate::types::{VariantSite, WindowId};

pub fn load_variants_json(path: &str, window: &WindowId) -> Result<Vec<VariantSite>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read variants file: {path}"))?;
    let all: Vec<VariantSite> = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse variants JSON: {path}"))?;

    Ok(all
        .into_iter()
        .filter(|v| v.pos >= window.start && v.pos <= window.end)
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
    fn loads_and_filters_variants_by_window() {
        let path = temp_path("variants");
        let json = r#"
[
  {"pos":101,"ref_allele":"A","alt_allele":"G","genotype":[0,1]},
  {"pos":251,"ref_allele":"C","alt_allele":"T","genotype":[0,1]},
  {"pos":9999,"ref_allele":"G","alt_allele":"A","genotype":[0,1]}
]
"#;
        std::fs::write(&path, json).unwrap();

        let window = WindowId {
            chrom: "chr20".to_string(),
            start: 100,
            end: 300,
        };

        let variants = load_variants_json(path.to_str().unwrap(), &window).unwrap();
        assert_eq!(variants.len(), 2);
        assert_eq!(variants[0].pos, 101);
        assert_eq!(variants[1].pos, 251);

        let _ = std::fs::remove_file(path);
    }
}
