nextflow.enable.dsl = 2

data = channel
        .fromPath(params.sample_alignments_tsv, type: "file", checkIfExists: true)
        .splitCsv(sep: "\t", header: ["sample_id", "alignment_file"])
        .map { row -> tuple(row.sample_id, row.alignment_file) }

reference_fasta = file(params.reference_fasta, type: "file", checkIfExists: true)
reference_fasta_fai = file("${reference_fasta}.fai", checkIfExists: true)

delly_exclude_regions_bed = file(params.delly_exclude_regions_bed, type: "file", checkIfExists: true)
smoove_exclude_regions_bed = file(params.smoove_exclude_regions_bed, type: "file", checkIfExists: true)
expansion_hunter_variant_catalog_json = file(params.expansion_hunter_variant_catalog_json, type: "file", checkIfExists: true)
recode_delly_python_script = file(params.recode_delly_python_script, type: "file", checkIfExists: true)
mosdepth_segmental_duplications_bed = file(params.mosdepth_segmental_duplications_bed, type:"file", checkIfExists: true)

// Run mosdepth for getting coverage statistics.
process COVERAGE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file) 
    output:
    path("${sample_id}.mosdepth.summary.txt")
    path("${sample_id}.segdups.regions.bed.gz") 
    script:
    """
    mosdepth --threads $task.cpus --no-per-base --fast-mode --fasta $reference_fasta $sample_id $alignment_file
    mosdepth --threads $task.cpus --by $mosdepth_segmental_duplications_bed --fasta $reference_fasta ${sample_id}.segdups $alignment_file
    """
    stub:
    """
    touch ${sample_id}.mosdepth.summary.txt
    """
}

// Run the Delly variant caller.
process DELLY {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-delly.bcf")
    script:
    """
    delly call -g $reference_fasta -x $delly_exclude_regions_bed  -o tmp.bcf $alignment_file -q 20 -s 15 -z 5
    delly call -g $reference_fasta -v tmp.bcf -x $delly_exclude_regions_bed -o $sample_id-delly.bcf  $alignment_file -q 20 -s 15 -z 5
    rm tmp.bcf
    """
    stub:
    """
    touch $sample_id-delly.bcf
    """
}
// Recode the Delly .bcf file and convert it into a .vcf file.
process RECODE {
    input:
    tuple val(sample_id), path(output_delly)
    output:
    tuple val(sample_id), path("$sample_id-delly-recode.vcf")
    script:
    """
    python $recode_delly_python_script $output_delly ${sample_id}-delly-recode.vcf
    """
    stub:
    """
    touch ${sample_id}-delly-recode.vcf
    """
}


// Run the Manta variant caller.
process MANTA {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-manta.vcf.gz")
    script:
    """
    configManta.py --bam $alignment_file --referenceFasta $reference_fasta --runDir run_folder/ 
    cd run_folder
    python runWorkflow.py
    mv results/variants/diploidSV.vcf.gz ../$sample_id-manta.vcf.gz
    """
    stub:
    """
    touch $sample_id-diploidSV.vcf.gz
    """
}

// Run the Smoove variant caller.
process SMOOVE {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), val(alignment_file)
    output:
    tuple val(sample_id), path("$sample_id-smoove.vcf.gz")
    script:
    """
    smoove call --name $sample_id --exclude $smoove_exclude_regions_bed --fasta $reference_fasta -p $task.cpus $alignment_file 
    """
    stub:
    """
    touch $sample_id-smoove.vcf.gz
    """
}

// Filter VCF with bcftools (filter on SVLEN and PASS flag).
process FILTER { 
    input:
    tuple val(sample_id), path(vcf)
    output: 
    tuple val(sample_id), path("${sample_id}_$vcf-filtered.vcf")
    script:
    if (vcf.name =~ /smoove/ )
    """                                                                         
    bcftools view --threads $task.cpus --include "SVLEN>=50 || SVLEN<=-50" $vcf > ${sample_id}_$vcf-filtered.vcf
    """
    else
    """
    bcftools view --threads $task.cpus --include "SVLEN>=50 || SVLEN<=-50" $vcf | bcftools view --threads $task.cpus --include "%FILTER='PASS'" > ${sample_id}_$vcf-filtered.vcf
    """
    stub:
    """
    touch ${sample_id}_$vcf-filtered.vcf
    """
}

// Run SURVIVOR to merge the within-sample variant calls.
process MERGE {
    input:
    val(sample_id)
    tuple path(vcf1), path(vcf2), path(vcf3)
    output: 
    tuple val(sample_id), path("${sample_id}_merged.vcf")
    script:
    """
    ls $vcf1    >  vcf_list.txt
    ls $vcf2    >> vcf_list.txt
    ls $vcf3    >> vcf_list.txt
    sort vcf_list.txt 
    SURVIVOR merge vcf_list.txt 1000 0 0 0 0 50 ${sample_id}_merged.vcf
    """
    stub:
    """
    echo $vcf1 >  vcf_list.txt
    echo $vcf2 >> vcf_list.txt
    echo $vcf3 >> vcf_list.txt
    sort vcf_list.txt
    touch ${sample_id}_merged.vcf
    """
}

process SORT {
    publishDir "results/$sample_id/", mode: "copy"
    input:
    tuple val(sample_id), path(vcf)
    output:
    path("${sample_id}.merged.sorted.vcf")
    script:
    """
    bcftools sort $vcf > ${sample_id}.merged.sorted.vcf
    """
}

 process EXPANSION_HUNTER {
    publishDir "results/$sample_id/", mode: "copy"
    input:                                                                   
    tuple val(sample_id), val(alignment_file)
    output: 
    path("${sample_id}.expansion_hunter.vcf")
    script:
    """
    ExpansionHunter --reads $alignment_file --reference $reference_fasta --variant-catalog $expansion_hunter_variant_catalog_json  --output-prefix ${sample_id}.expansion_hunter --analysis-mode streaming --threads $task.cpus
    """
    stub:
    """
    touch ${sample_id}.eh.vcf
    """
 }


workflow {
    // Get coverage using Mosdepth 
    COVERAGE(data)

    // Run variant calling (using Delly, Manta, and Smoove)
    delly_output = RECODE(DELLY(data))
    manta_output = MANTA(data)
    smoove_output = SMOOVE(data)

    // Filter all the VCFs, and then group tuples (by sample ID), then split into 2 channels using multimap (the first channel is the sample ID, and the second channel is the tuple of filtered VCFs for each sample)
    filtered_vcf_output = FILTER(delly_output.mix(manta_output, smoove_output))
    filtered_vcf_output.groupTuple().multiMap{ it -> sample_id: it[0]; vcfs: it[1] }.set{ filtered_vcfs_to_merge }

    // Merge the VCFs sample-wise using SURVIVOR (and then sort the output VCFs)
    merged_vcf_output = SORT(MERGE(filtered_vcfs_to_merge.sample_id, filtered_vcfs_to_merge.vcfs))

     expansion_hunter_output = EXPANSION_HUNTER(data)
    // hipstr_output = HIPSTR(data)
}
