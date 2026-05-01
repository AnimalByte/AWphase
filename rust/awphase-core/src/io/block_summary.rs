use std::fs;

use anyhow::{Context, Result};

use crate::bridge::block_summary::BlockSummary;

pub fn load_block_summaries_json(path: &str) -> Result<Vec<BlockSummary>> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read block summary file: {path}"))?;
    let blocks: Vec<BlockSummary> = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse block summary JSON: {path}"))?;
    Ok(blocks)
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
    fn loads_block_summaries() {
        let path = temp_path("blocks");
        let json = r#"
[
  {
    "block_id": "b1",
    "orientation_bias": 0.9,
    "support_weight": 3.2,
    "informative_sites": 4,
    "abstained_sites": 0,
    "mean_confidence": 0.8,
    "signed_confidence_sum": 3.2,
    "coherence": 1.0
  },
  {
    "block_id": "b2",
    "orientation_bias": -0.8,
    "support_weight": 3.0,
    "informative_sites": 4,
    "abstained_sites": 0,
    "mean_confidence": 0.75,
    "signed_confidence_sum": -3.0,
    "coherence": 1.0
  }
]
"#;
        std::fs::write(&path, json).unwrap();

        let blocks = load_block_summaries_json(path.to_str().unwrap()).unwrap();
        assert_eq!(blocks.len(), 2);
        assert_eq!(blocks[0].block_id, "b1");
        assert_eq!(blocks[1].block_id, "b2");

        let _ = std::fs::remove_file(path);
    }
}
