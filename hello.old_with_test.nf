#!/usr/bin/env nextflow     
nextflow.enable.dsl = 1

//params.alignment_file = "/expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/REACH000236.chr22.cram"
//params.alignment_file_index = "/expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/REACH000236.chr22.cram.crai"
// reference_file = "/expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/Homo_sapiens_assembly38.fasta"
// reference_file_index = "/expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/Homo_sapiens_assembly38.fasta.fai"
//params.reference_file = "/expanse/lustre/projects/ddp195/eiovino/pipeline_singularity_bind_path_tests/Homo_sapiens_assembly38.fasta"
//params.reference_file_index = "/expanse/lustre/projects/ddp195/eiovino/pipeline_singularity_bind_path_tests/Homo_sapiens_assembly38.fasta.fai"


//alignment_file = params.alignment_file
//alignment_file_index = params.alignment_file_index
//reference_file = params.reference_file
//reference_file_index = params.reference_file_index

process mosdepth_test {
     output: 
     file "output/version.txt" into mosdepth_file
     //file "${ccs_bam.simpleName}.ccs.fastq" into ccs_fastq_file
     script:
     """
     mkdir output
     mosdepth --version > output/version.txt
     """
}

process bcftools_test {                                                       
     script:                                                                  
     """                                                                                                                                   
     bcftools --version > version.txt                                  
     """                                                                      
}  

process output_test {
    input:
    file mosdepth_output from mosdepth_file

    script:
    """
    echo $mosdepth_output > my_version.txt
    """
}

 Delly singularity container does not have bash...
 process delly_test {
     script:
     """
     delly call -g $reference_file -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
     """
 }

 process expansionhunter_test {                                                         
      script:                                                                  
      """                                                                      
      ExpansionHunter --help
      """                                                                      
 }           

process survivor_test {                                                         
     script:                                                                  
     """                                                                      
     SURVIVOR &> ls.txt                                                             
     """                                                                      
}  

process smoove_test {
    input:
    alignment_file
    alignment_file_index
    reference_file
    reference_file_index

    script:
    """
    smoove call --name REACH000236 --fasta $reference_file $alignment_file
    """
}

process manta_test {
     output: 
     
     script:
     """
     configManta.py --bam $alignment_file --referenceFasta $reference_file --runDir run_folder/
     # cd run_folder
     # python runWorkflow.py
     """
}

  params.input_csv = "input_data.csv"
  Channel
       .fromPath(params.input_csv)
       .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
       .map { row -> tuple(row.sample_id, row.alignment_file) }
       .into {manta_channel ; delly_channel ; ExpansionHunter_channel ;smoove_channel}
//   
//  // // // alignment_file = file("/expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/REACH000236.chr22.cram")
  reference_file = file("/expanse/lustre/projects/ddp195/eiovino/fasta/Homo_sapiens_assembly38.fasta")
  // 
  process runManta {
       input:
       set sample_id, alignment_file from manta_channel
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
    input:                                                                      
    file output_manta from manta_file                                     
    file output_index_manta from manta_index_file
    output:
    file "diploidSV_filtered.vcf" into manta_output_file                                                                            
    script:                                                                     
    """                                                                         
    bcftools view -i 'SVLEN>=50 || SVLEN<=-50' $output_manta > tmp.vcf
    bcftools view -i "%FILTER='PASS'" tmp.vcf > diploidSV_filtered.vcf
    """
}  

  process runDelly {
       input: 
       set sample_id, alignment_file from delly_channel
       output:
       file "$sample_id-tmp-delly.vcf" into delly_file 
       script:
        """
        delly call -g $reference_file -x /expanse/lustre/projects/ddp195/j3guevar/sv_testing/resources/human.hg38.excl.tsv  -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
        delly call -g $reference_file -v tmp.bcf -x /expanse/lustre/projects/ddp195/j3guevar/sv_testing/resources/human.hg38.excl.tsv -o $sample_id-delly.bcf  $alignment_file -q 20 -s 15 -z 5
        bcftools view $sample_id-delly.bcf > $sample_id-tmp-delly.vcf
        rm $sample_id-delly.bcf tmp.bcf
        """
}

process post_processingDellyTest {                                                            
       input:  
       file output_delly from delly_file
       output:
       file "final-delly.vcf" into delly_output_file
       script:
       """
       python /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/procesing_delly.py $output_delly
       bcftools view -i 'SVLEN>=50 ' out.vcf > tmpDelly.vcf
       bcftools view -i "%FILTER='PASS'" tmpDelly.vcf > final-delly.vcf
       rm tmpDelly.vcf out.vcf
       """
}


  
  process runSmoove {                                                          
       input:                                                                     
       set sample_id, alignment_file from smoove_channel
       output:
       file "$sample_id-smoove.vcf.gz" into smoove_file                           
       script:                                                                    
        """                                                                       
        smoove call --name $sample_id --fasta $reference_file -p 8 $alignment_file 
        """                                                                       
  //                                                                              
   } 

process post_processing_smmove {                                                           
    input:                                                                      
    file output_smoove from smoove_file                                           
    output: 
    file "filtered_smoove.vcf" into smoove_ouput_file
    script:                                                                     
    """                                                                         
    bcftools view -i 'SVLEN>=50 || SVLEN<=-50' $output_smoove > filtered_smoove.vcf          
    """                                                                         
}  

process runExpansionHunter {                                                           
       input:                                                                   
       set sample_id, alignment_file from ExpansionHunter_channel                                                 
       script:                                                                  
       """                                                                       
       ExpansionHunter --reads $alignment_file --reference $reference_file --variant-catalog /expanse/lustre/projects/ddp195/eiovino/sv_pipeline_nf/data/variant_catalog.json  --output-prefix $sample_id
       """                                                                                                                                                 
   } 


process runSurvivor {                                                    
       input:
       file filtered_manta from manta_output_file                                                                    
       file filtered_smoove from smoove_ouput_file
       file filtered_delly from delly_output_file              
       script:                                                                  
       """                                                                       
       ls $filtered_manta   > vcf_list.txt 
       ls $filtered_smoove  >> vcf_list.txt
       ls $filtered_delly   >> vcf_list.txt
       SURVIVOR merge vcf_list.txt 1000 1 1 1 1 50 sample_merged.vcf 
       """                                                                                                                         
   }                                                                            
         


 
