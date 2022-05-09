#!/usr/bin/env nextflow     
nextflow.enable.dsl = 1
params.help = false
if (params.help) {
    log.info """
    -----------------------------------------------------------------------
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

params.excluded_regions_delly = "Hg38/exlude.regions.delly.human.hg38.excl.tsv"
excluded_regions_delly = file(params.excluded_regions_delly, type: "file") 

outdir= params.outdir ?: './results'

  Channel
       .fromPath(alignment_list)
       .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
       .map { row -> tuple(row.sample_id, row.alignment_file) }
       .into {manta_channel ; manta_post_channel ;delly_channel ; delly_post_channel; ExpansionHunter_channel ; 
       smoove_channel; smoove_post_channel; mosdepth_channel}
  
  process runManta {
       input:
       set sample_id, alignment_file from manta_channel
       file reference_file_index
       output:
       file "run_folder/results/variants/diploidSV.vcf.gz" into manta_file
       file "run_folder/results/variants/diploidSV.vcf.gz.tbi" into manta_index_file
  
       script:
       """
       configManta.py --bam $alignment_file --referenceFasta $reference_file --runDir run_folder/ 
       cd run_folder
       python runWorkflow.py
       """
   }

process manta_post_processing {                                                           
    publishDir "$outdir/results/manta_out", mode: 'copy'

    input:
    set sample_id, alignment_file from manta_post_channel                                                                      
    file output_manta from manta_file                                     
    file output_index_manta from manta_index_file
    output:
    file "$sample_id-manta_filtered.vcf" into manta_output_file                                                                            
    script:                                                                     
    """                                                                         
    bcftools view -i 'SVLEN>=50 || SVLEN<=-50' $output_manta > tmp.vcf
    bcftools view -i "%FILTER='PASS'" tmp.vcf > $sample_id-manta_filtered.vcf
    """
}  

  process runDelly {
       input: 
       set sample_id, alignment_file from delly_channel
       file reference_file_index
       output:
       file "$sample_id-tmp-delly.vcf" into delly_file 
       script:
        """
        delly call -g $reference_file -x $excluded_regions_delly  -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
        delly call -g $reference_file -v tmp.bcf -x $excluded_regions_delly -o $sample_id-delly.bcf  $alignment_file -q 20 -s 15 -z 5
        bcftools view $sample_id-delly.bcf > $sample_id-tmp-delly.vcf
        rm $sample_id-delly.bcf tmp.bcf
        """
}

process post_processingDellyTest {
       publishDir "$outdir/results/delly_out", mode: 'copy' 
       input:
       set sample_id, alignment_file from delly_post_channel  
       file output_delly from delly_file
       output:
       file "$sample_id-final-delly.vcf" into delly_output_file
       script:
       """
       python /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/procesing_delly.py $output_delly
       bcftools view -i 'SVLEN>=50 ' out.vcf > tmpDelly.vcf
       bcftools view -i "%FILTER='PASS'" tmpDelly.vcf > $sample_id-final-delly.vcf
       rm tmpDelly.vcf out.vcf
       """
}


  
  process runSmoove {                                                          
       input:                                                                     
       set sample_id, alignment_file from smoove_channel
       file reference_file_index
       output:
       file "$sample_id-smoove.vcf.gz" into smoove_file                           
       script:                                                                    
        """                                                                       
        smoove call --name $sample_id --fasta $reference_file -p 8 $alignment_file 
        """                                                                       
                                                                              
   } 

process post_processing_smoove { 
    publishDir "$outdir/results/smoove_out", mode: 'copy' 
    input:
    set sample_id from smoove_post_channel                                                                      
    file output_smoove from smoove_file                                           
    output: 
    file "$sample_id-filtered_smoove.vcf" into smoove_ouput_file
    script:                                                                     
    """                                                                         
    bcftools view -i 'SVLEN>=50 || SVLEN<=-50' $output_smoove > $sample_id-filtered_smoove.vcf          
    """                                                                         
}  

process runExpansionHunter {
       publishDir "$outdir/results/ExpansionHunter_out", mode: 'copy'
       input:                                                                   
       set sample_id, alignment_file from ExpansionHunter_channel
       output: 
       file "*.vcf" into ExH_output
       script:                                                                  
       """                                                                       
       ExpansionHunter --reads $alignment_file --reference $reference_file --variant-catalog /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/Hg38/ExpansioHunter_variant_catalog.json  --output-prefix $sample_id
       """                                                                                                                                                 
   } 


process runSurvivor {                                                    
       publishDir "$outdir/results/survivor_out", mode: 'copy'
       input:
       file filtered_manta from manta_output_file                                                                    
       file filtered_smoove from smoove_ouput_file
       file filtered_delly from delly_output_file              
       output: 
       file "sample_merged.vcf" into survivor_output 
       script:                                                                  
       """                                                                       
       ls $filtered_manta   > vcf_list.txt 
       ls $filtered_smoove  >> vcf_list.txt
       ls $filtered_delly   >> vcf_list.txt
       SURVIVOR merge vcf_list.txt 1000 1 1 1 1 50 sample_merged.vcf 
       """                                                                                                                         
   }                                                                            
         
//Mosdepth running NOTE: Add script for plotting te cov result? 
process mosdepth {                                                           
       publishDir "$outdir/results/mosdepth_out", mode: 'copy'
       input:                                                                     
       set sample_id, alignment_file from mosdepth_channel
       output: 
       file "*mosedepth.wgs.mosdepth.*" into mosdepth_out                                                                    
                                                                                
       script:                                                                  
       """                                                                       
       mosdepth -n --fast-mode --by 500 $sample_id-mosedepth.wgs  $alignment_file --fasta $reference_file
       """                                                                      
   } 

 
