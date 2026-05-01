use std::fs;
use std::io::{Cursor, Write};

use anyhow::{Context, Result};

use crate::candidate::interface::PanelHaplotype;

pub fn load_panel_haplotypes(path: &str) -> Result<Vec<PanelHaplotype>> {
    if path.ends_with(".bin") {
        let bytes = fs::read(path)
            .with_context(|| format!("failed to read panel haplotype binary file: {path}"))?;
        let mut cur = Cursor::new(bytes);
        let panel: Vec<PanelHaplotype> = bincode::deserialize_from(&mut cur)
            .with_context(|| format!("failed to parse panel haplotype binary: {path}"))?;
        Ok(panel)
    } else {
        load_panel_haplotypes_json(path)
    }
}

pub fn load_panel_haplotypes_json(path: &str) -> Result<Vec<PanelHaplotype>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read panel haplotype file: {path}"))?;
    let panel: Vec<PanelHaplotype> = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse panel haplotype JSON: {path}"))?;
    Ok(panel)
}

pub fn save_panel_haplotypes_bin(path: &str, panel: &[PanelHaplotype]) -> Result<()> {
    let bytes = bincode::serialize(panel)
        .with_context(|| format!("failed to serialize panel haplotypes for {path}"))?;
    let mut fh = fs::File::create(path)
        .with_context(|| format!("failed to create panel haplotype binary file: {path}"))?;
    fh.write_all(&bytes)
        .with_context(|| format!("failed to write panel haplotype binary file: {path}"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use std::time::{SystemTime, UNIX_EPOCH};

    fn temp_path(name: &str, ext: &str) -> PathBuf {
        let ts = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        std::env::temp_dir().join(format!("awphase_{name}_{ts}.{ext}"))
    }

    #[test]
    fn loads_panel_haplotypes_json_roundtrip() {
        let path = temp_path("panel_json", "json");
        let json = r#"
[
  {
    "donor_id": "d1",
    "hap_id": 0,
    "group": "eur",
    "alleles": [1,1,-1,-1]
  },
  {
    "donor_id": "d2",
    "hap_id": 1,
    "group": "afr",
    "alleles": [1,-1,1,-1]
  }
]
"#;
        std::fs::write(&path, json).unwrap();

        let panel = load_panel_haplotypes_json(path.to_str().unwrap()).unwrap();
        assert_eq!(panel.len(), 2);
        assert_eq!(panel[0].donor_id, "d1");
        assert_eq!(panel[1].group, "afr");

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn loads_panel_haplotypes_binary_roundtrip() {
        let path = temp_path("panel_bin", "bin");
        let panel = vec![
            PanelHaplotype {
                donor_id: "d1".to_string(),
                hap_id: 0,
                group: "eur".to_string(),
                alleles: vec![1, 1, -1, -1],
            },
            PanelHaplotype {
                donor_id: "d2".to_string(),
                hap_id: 1,
                group: "afr".to_string(),
                alleles: vec![1, -1, 1, -1],
            },
        ];

        save_panel_haplotypes_bin(path.to_str().unwrap(), &panel).unwrap();
        let loaded = load_panel_haplotypes(path.to_str().unwrap()).unwrap();

        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded[0].donor_id, "d1");
        assert_eq!(loaded[1].group, "afr");

        let _ = std::fs::remove_file(path);
    }
}
