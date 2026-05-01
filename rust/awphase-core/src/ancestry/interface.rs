use anyhow::Result;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroupWeight {
    pub label: String,
    pub weight: f32,
}

pub trait AncestryModel {
    fn get_global_prior(&self, sample_id: &str) -> Result<Vec<GroupWeight>>;
}
