#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    selectionsites
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/vmiro/Selection_sites
    Author : Veronica Miro Pina and Claudia Serrano Colome
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// SUBWORKFLOW: Consisting of regressions and simulations per bin type
//

include { REG_SIM } from '/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/subworkflows/reg_sim'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



cancertypes = Channel.fromPath('cancers_3mers.txt').splitText().map{ it -> it.trim() }
//effect_considered = Channel.from("all", "synonymous", "nonsynonymous")  
effect_considered = Channel.from("all") 
cancer_effect = cancertypes.combine(effect_considered)

process WAIT{
    label 'bashprocess'

    input:
    val cancer from cancer_effect_in.collect()

    output:
    val cancer

    script:
    """
    echo ${cancer}
    """
   }
   
process COMPARISONS{
    label 'pythonprocess'

    input:
    val cancer_effect_in

    output:
    val cancer_effect_in

    script:
    """
    python ${baseDir}/bin/comparison.py ${cancer_effect_in[0]} ${cancer_effect_in[1]} 
    """
   }

process FINAL{
    label 'bashprocess'

    input:
    val cancer_effect_in

    output:
    stdout

    script:
    """
    cat nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/best/chosen_model_${cancer_effect_in[0]}.csv
    """
   }
 

process GZIP{
    label 'bashprocess'

    input:
    val cancer_effect_in

    output:
    stdout

    script:
    """
    gzip nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/${cancer_effect_in[0]}/simulated_data/*
    """
   }


workflow {
     REG_SIM(cancer_effect)
     wait_out = WAIT(REG_SIM.out)
     com_out = COMPARISONS(wait_out)
     FINAL(com_out)
     GZIP(com_out)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
