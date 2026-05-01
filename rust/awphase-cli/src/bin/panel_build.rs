use anyhow::Result;
use awphase_core::io::panel::{load_panel_haplotypes_json, save_panel_haplotypes_bin};
use clap::Parser;

#[derive(Parser, Debug)]
struct Cli {
    #[arg(long)]
    input: String,
    #[arg(long)]
    output: String,
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    let panel = load_panel_haplotypes_json(&cli.input)?;
    save_panel_haplotypes_bin(&cli.output, &panel)?;
    eprintln!("wrote {} haplotypes to {}", panel.len(), cli.output);
    Ok(())
}
