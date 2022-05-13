#!/usr/bin/env nextflow     
nextflow.enable.dsl = 1
params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------
    Required arguments: 
    --reference_file    Absolute path of reference FASTA file. Index file must be in the same directory
    --alignment_list    Tab delimited file with two column: sample_id and absolute path of aligned files. The index file must be same directory
    --outdir            Directory for output file
    """
}


reference_file = file(params.reference_file, type: "file")
if (!reference_file.exists()) {
    println("--reference_file: File doesn't exist, check path.")
    exit 1
}

reference_file_index = file("${reference_file}.fai")
if (!reference_file_index.exists()) {
    println("Missing Fasta Index")                
    exit 1                                                      
}


alignment_list = file(params.alignment_list, type:"file")
if (!alignment_list.exists()){
   println ("--alignment_list: File doesn't exist, check alignment file list.")
   exit 1
}

project_dir = projectDir
params.excluded_regions_delly = "$project_dir/Hg38/exlude.regions.delly.human.hg38.excl.tsv"
excluded_regions_delly = file(params.excluded_regions_delly, type: "file") 
params.excluded_regions_smoove = "$project_dir/Hg38/exclude.smoove.cnvnator_100bp.GRCh38.20170403.bed" 
excluded_regions_smoove = file(params.excluded_regions_smoove, type: "file")  
params.variant_catalog = "$project_dir/Hg38/ExpansioHunter_variant_catalog.json"
variant_catalog = file(params.variant_catalog, type: "file") 
println "Project : $workflow.projectDir"
outdir= params.outdir ?: './results'
log.info("\n")
log.info(println "Project : $workflow.projectDir" )
log.info("\n")

Channel
     .fromPath(alignment_list)
     .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
     .map { row -> tuple(row.sample_id, row.alignment_file) }
     .into {manta_channel; delly_channel; ExpansionHunter_channel; smoove_channel; mosdepth_channel}
  
//println "Sample : $sample_id , $alignment_file" 

process runManta {
     input:
     tuple val(sample_id), val(alignment_file) from manta_channel
     // set sample_id, alignment_file from manta_channel
     file reference_file_index
     output:
     tuple val(sample_id), file("run_folder/results/variants/$sample_id-diploidSV.vcf.gz"), file("run_folder/results/variants/$sample_id-diploidSV.vcf.gz.tbi") into manta_output_channel
     // file "run_folder/results/variants/$sample_id-diploidSV.vcf.gz" into manta_file
     // file "run_folder/results/variants/$sample_id-diploidSV.vcf.gz.tbi" into manta_index_file
 
     script:
     """
     configManta.py --bam $alignment_file --referenceFasta $reference_file --runDir run_folder/ 
     cd run_folder
     python runWorkflow.py
     mv results/variants/diploidSV.vcf.gz results/variants/$sample_id-diploidSV.vcf.gz
     mv results/variants/diploidSV.vcf.gz.tbi results/variants/$sample_id-diploidSV.vcf.gz.tbi
     """
 }

 process manta_post_processing {                                                           
     publishDir "$outdir/$sample_id", mode: 'copy'
 
     input:
     tuple val(sample_id), file(output_manta), file(output_index_manta) from manta_output_channel
     // file output_manta from manta_file                                     
     // file output_index_manta from manta_index_file
     output:
     // file "$output_manta-manta_filtered.vcf" into manta_output_file                                                                            
     tuple val(sample_id), file("$output_manta-manta_filtered.vcf") into manta_output_channel2
     script:                                                                     
     """                                                                         
     bcftools view -i 'SVLEN>=50 || SVLEN<=-50' $output_manta > $output_manta-tmp.vcf
     bcftools view -i "%FILTER='PASS'" $output_manta-tmp.vcf > $output_manta-manta_filtered.vcf
     """
 }  

process runDelly {
    input: 
    //set sample_id, alignment_file from delly_channel
    tuple val(sample_id), val(alignment_file) from delly_channel
    file reference_file_index
    
    output:
    // file "$sample_id-delly.bcf" into delly_file 
    tuple val(sample_id), file("$sample_id-delly.bcf") into delly_output_channel

    script:
    """
    delly call -g $reference_file -x $excluded_regions_delly  -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
    delly call -g $reference_file -v tmp.bcf -x $excluded_regions_delly -o $sample_id-delly.bcf  $alignment_file -q 20 -s 15 -z 5
    rm tmp.bcf 
    """
}

process python_delly {
    input:
    tuple val(sample_id), file(output_delly) from delly_output_channel
    // file output_delly  from delly_file
    output:
    // file "$output_delly-recode.vcf"  into python_output
    tuple val(sample_id), file("$output_delly-recode.vcf") into delly_output_channel2
    script:
    """
    python $project_dir/processing_delly.py $output_delly $output_delly-recode.vcf
    """
}

process post_processing_delly {
    publishDir "$outdir/$sample_id", mode: 'copy'
    input:                                                                   
    // file output_python from python_output
    tuple val(sample_id), file(output_python) from delly_output_channel2
    output:
    // file "$output_python-final_delly.vcf" into delly_output_file
    tuple val(sample_id), file("$output_python-final_delly.vcf") into delly_output_channel3
    script:
    """
    bcftools view -i 'SVLEN>=50 ' $output_python | bcftools view -i "%FILTER='PASS'" > $output_python-final_delly.vcf
    """
}


  
process runSmoove {                                                          
     input:                                                                     
     // set sample_id, alignment_file from smoove_channel
     tuple val(sample_id), val(alignment_file) from smoove_channel
     file reference_file_index
     output:
     // file "$sample_id-smoove.vcf.gz" into smoove_file                           
     tuple val(sample_id), file("$sample_id-smoove.vcf.gz") into smoove_output_channel
     script:                                                                    
     """                                                                       
     smoove call --name $sample_id --exclude $excluded_regions_smoove --fasta $reference_file -p 8 $alignment_file 
     """                                                                       
} 

process post_processing_smoove { 
    publishDir "$outdir/$sample_id", mode: 'copy' 
    input:
    tuple val(sample_id), file(output_smoove) from smoove_output_channel 
    // file output_smoove from smoove_file                                           
    output: 
    // file "$output_smoove-filtered.vcf" into smoove_output_file
    tuple val(sample_id), file("$output_smoove-filtered.vcf") into smoove_output_channel2
    script:                                                                     
    """                                                                         
    bcftools view -i 'SVLEN>=50 || SVLEN<=-50' $output_smoove > $output_smoove-filtered.vcf          
    """                                                                         
}  

//process runExpansionHunter {
//       publishDir "$outdir/$sample_id", mode: 'copy'
//       input:                                                                   
//       set sample_id, alignment_file from ExpansionHunter_channel
//       output: 
//       file "*.vcf" into ExH_output
//       script:                                                                  
//       """                                                                       
//       ExpansionHunter --reads $alignment_file --reference $reference_file --variant-catalog $variant_catalog  --output-prefix $sample_id
//       """                                                                                                                                                 
//   } 


process runSurvivor {                                                    
       publishDir "$outdir/$sample_id", mode: 'copy'
       input:
       tuple val(sample_id), file(filtered_manta), file(filtered_smoove), file(filtered_delly) from manta_output_channel2.join(smoove_output_channel2).join(delly_output_channel3)
       //file filtered_manta from manta_output_file                                                                    
       //file filtered_smoove from smoove_output_file
       //file filtered_delly from delly_output_file             
       output: 
       file "${sample_id}_merged.test.vcf" into survivor_output 
       script:                                                                  
       """                                                                       
       ls $filtered_manta   > vcf_list.txt 
       ls $filtered_smoove  >> vcf_list.txt
       ls $filtered_delly   >> vcf_list.txt
       SURVIVOR merge vcf_list.txt 1000 1 1 1 1 50 ${sample_id}_merged.test.vcf 
       """                                                                                                                         
}                                                                            
         
 process mosdepth {                                                           
        publishDir "$outdir/$sample_id", mode: 'copy'
        input:                                                                     
        set sample_id, alignment_file from mosdepth_channel
        file reference_file_index
        output: 
        file "*mosedepth.wgs.mosdepth.*" into mosdepth_out                                                                    
                                                                                 
        script:                                                                  
        """                                                                       
        mosdepth -n --fast-mode --by 500 $sample_id-mosedepth.wgs -f $reference_file $alignment_file
        """                                                                      
 } 

 
