nextflow.enable.dsl = 2

params.label = params.label ?: 'chr22_20_25mb'

workflow {
    throw new IllegalStateException(
        "AWPhase does not yet have a production Nextflow pipeline. " +
        "Use `pixi run phase9d-window plan ${params.label}` to inspect the " +
        "current manifest-driven Phase9D path without downloading data."
    )
}
