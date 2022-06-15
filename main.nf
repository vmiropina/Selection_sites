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
effect_considered = Channel.from("all", "synonymous", "nonsynonymous")  
cancer_effect = cancertypes.combine(effect_considered)
   
process COMPARISONS{
    label 'pythonprocess1'

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
    var="\$(cat /nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/best/chosen_model_align_${cancer_effect_in[0]}.csv)"
    cp /nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/${cancer_effect_in[0]}/plots/hist_\${var}.png /nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/best/hist_align_${cancer_effect_in[0]}.png
    var="\$(cat /nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/best/chosen_model_${cancer_effect_in[0]}.csv)"
    cp /nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/${cancer_effect_in[0]}/plots/hist_\${var}.png /nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/best/hist_noalign_${cancer_effect_in[0]}.png
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
    gzip /nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/syn/${cancer_effect_in[1]}/${cancer_effect_in[0]}/simulated_data/*
    """
   }


workflow {
     REG_SIM(cancer_effect)
     com_out = COMPARISONS(REG_SIM.out.collect())
     FINAL(com_out)
     //GZIP(com_out)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
