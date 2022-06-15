//
// Do regressions and simulations per bintype
//



bintypes = Channel.fromPath('bins.txt').splitText().map{ it -> it.trim() }

process REGRESSIONS{
    label 'rprocess'

    input:
    val cancertypes_in

    output:
    val cancertypes_in

    script:
    """
    Rscript ${baseDir}/bin/regressions.R ${cancertypes_in[0]} ${cancertypes_in[2]} ${cancertypes_in[1]}
    """
   }

process SIMULATIONS{
    label 'pythonprocess1'

    input:
    val cancertypes_in
    
    output:
    val cancertypes_in

    script:
    """
    python ${baseDir}/bin/simulate_selection_sites.py ${cancertypes_in[0]} ${cancertypes_in[2]} ${cancertypes_in[1]}
    """
   }

process PLOTS{
    label 'pythonprocess2'

    input:
    val cancertypes_in
    
    output:
    val cancertypes_in

    script:
    """
    python ${baseDir}/bin/plots.py ${cancertypes_in[0]} ${cancertypes_in[2]} ${cancertypes_in[1]}
    """
   }
   




workflow REG_SIM {
    take:
    cancertype // cancertype read in the main workflow
    
    main:
    reg_out = REGRESSIONS(cancertype.combine(bintypes))
    simu_out = SIMULATIONS(reg_out)
    plot_out = PLOTS(simu_out)

    emit:
    plot_out
}

