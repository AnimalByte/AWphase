use std::fs;

use anyhow::{Context, Result};

use crate::types::{SiteDonorPriors, WindowId};

pub fn load_site_priors_json(path: &str, window: &WindowId) -> Result<SiteDonorPriors> {
    let text = fs::read_to_string(path)
        .with_context(|| format!("failed to read site priors file: {path}"))?;
    let mut all: SiteDonorPriors = serde_json::from_str(&text)
        .with_context(|| format!("failed to parse site priors JSON: {path}"))?;

    all.entries
        .retain(|e| e.pos >= window.start && e.pos <= window.end);
    all.entries.sort_by_key(|e| e.pos);

    Ok(all)
}

pub fn site_bias_for_pos(priors: &SiteDonorPriors, pos: u64) -> f32 {
    match priors.entries.binary_search_by_key(&pos, |e| e.pos) {
        Ok(idx) => priors.entries[idx].bias,
        Err(_) => 0.0,
    }
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
    fn loads_and_filters_site_priors_by_window() {
        let path = temp_path("site_priors");
        let json = r#"
{
  "entries": [
    { "pos": 101, "bias": 0.2 },
    { "pos": 251, "bias": -0.3 },
    { "pos": 9999, "bias": 0.9 }
  ]
}
"#;
        std::fs::write(&path, json).unwrap();

        let window = WindowId {
            chrom: "chr20".to_string(),
            start: 100,
            end: 300,
        };

        let priors = load_site_priors_json(path.to_str().unwrap(), &window).unwrap();
        assert_eq!(priors.entries.len(), 2);
        assert_eq!(priors.entries[0].pos, 101);
        assert_eq!(priors.entries[0].bias, 0.2);
        assert_eq!(priors.entries[1].pos, 251);
        assert_eq!(priors.entries[1].bias, -0.3);

        let _ = std::fs::remove_file(path);
    }

    #[test]
    fn site_bias_lookup_returns_expected_values() {
        let priors = SiteDonorPriors {
            entries: vec![
                crate::types::SitePriorEntry {
                    pos: 101,
                    bias: 0.2,
                },
                crate::types::SitePriorEntry {
                    pos: 251,
                    bias: -0.3,
                },
            ],
        };

        assert_eq!(site_bias_for_pos(&priors, 101), 0.2);
        assert_eq!(site_bias_for_pos(&priors, 251), -0.3);
        assert_eq!(site_bias_for_pos(&priors, 999), 0.0);
    }
}
