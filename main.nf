params.help = false

if (params.help) {
    log.info """
    SV calling workflow for short-read sequencing data
    Required arguments:
    --reference
    --cram
    """.stripIndent()
    exit 0
}


process run_erds {
    script:
    """
    perl erds_pipeline.pl -r $reference -b $alignment_file -v $snv_vcf_file -o ${alignment_file.simpleName}
    """
}


// For Delly, we'll need to postprocess the output 
process run_delly {
    out:
    "${input_cram.simpleName}.bcf" into delly_bcf

    script:
    """
    delly call --exclude $delly_exclude_file --genome $reference --outfile ${alignment_file.simpleName}.bcf --map-qual 20 --mad-cutoff 15 --min-clique-size 5 $alignment_file
    """
}

process fix_delly_output {
    in:
    delly_bcf_input from delly_bcf

    script:
    """
    bcftools view -Oz -o delly_output.vcf.gz $delly_bcf_input
    tabix delly_output.vcf.gz
    """
} 


process run_smoove {
    out:

    script:
    """
    smoove call --outdir results-smoove/ --exclude $smoove_exclude_file --name ${alignment_file.simpleName} --fasta $reference --processes 1 --genotype $alignment_file
    """
}

process run_manta {
    out:
    "run_folder/results/variants/diploidSV.vcf.gz" as manta_vcf
    "run_folder/results/variants/diploidSV.vcf.gz.tbi" as manta_vcf_index

    script:
    """
    configManta.py --bam $alignment_file --referenceFasta $reference --runDir run_folder/ && python run_folder/runWorkflow.py
    """
}

// For SVaba, we'll need to postprocess the output
process run_svaba {

    script:
    """
    svaba run -t $alignment_file -p $threads -I -a ${alignment_file.simpleName} -G $reference
    """
}

process fix_svaba_output {
    script:
    """
    convert_SvABA_vcf.original.py $svaba_input_vcf > $svaba_output_vcf
    """
}


process run_expansion_hunter {
    script:
    """
    ExpansionHunter --reads $alignment_file --reference $reference --variant-catalog $expansion_hunter_variant_catalog_file --output-prefix ${alignment_file.simpleName}
    """
}

process run_hipstr {

    script:
    """
    HipSTR --chrom chr1 --min-reads 5 --def-stutter-model --bams $alignment_file --fasta $reference --regions $hipstr_bed_file --str-vcf ${alignment_file.simpleName}.hipstr.chr1.vcf.gz --bam-samps ${alignment_file.simpleName} --bam-libs ${alignment_file.simpleName}
    """
}

process run_cnvnator {
    script:
    """
    """
}

// Try 500bp distance rather than 1kb distance?
process run_survivor {
    script:
    """
    ls $erds_vcf_input    > vcf_list.txt
    ls $delly_vcf_input  >> vcf_list.txt
    ls $smoove_vcf_input >> vcf_list.txt
    ls $manta_vcf_input  >> vcf_list.txt
    ls $svaba_vcf_input  >> vcf_list.txt
    SURVIVOR merge vcf_list.txt 1000 1 1 1 1 50 sample_merged.vcf 
    """
}

process convert_survivor_vcf_to_bed {
    script:
    """
    bedtools sort -i $survivor_vcf_input > $survivor_bed_output
    """
}

process run_meta_classifier {
    script:
    """
    python run_meta_classifier.py --input_bed $survivor_bed_input --output_bed $classifier_bed_output
    """
}

process run_sv2_genotyper {
    script:
    """
    python run_sv2_genotyper.py --reference $reference --alignment_file $alignment_file_input --input_bed $classifier_bed_input --snv_vcf $snv_vcf_input 
    """
}
