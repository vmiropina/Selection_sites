process REGRESSIONS {
    tag "$regressions"


    input:
    path samplesheet

    output:
    path '*.csv'       , emit: csv
    path "versions.yml", emit: versions

    script: // This script is bundled with the pipeline, in nf-core/selectionsites/bin/
    """
    Rscript scripts/covariates_regressions/PCA_and_regression.R ${cancer_type} ${bin} ${effect}S
    """
}
