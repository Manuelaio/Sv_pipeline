import pysam

vcf_iterator = pysam.VariantFile("ita_cohort_survivor_intersect_merged_ACMG.vcf.gz")
vcf_out = pysam.VariantFile("ita_cohort_survivor_intersect_merged_ACMG_scored.vcf.gz", mode = "w", header = vcf_iterator.header)
vcf_out.header.add_meta(key = "FORMAT", items = [ ("ID", "SCORE"), ("Number", "."), ("Type", "Integer"), ("Description", "SCORE Classification using callers and length {ZERO: 0, WEAK: 1, MEDIUM: 2, STRONG: 3}")])

index_sample_dict = dict(enumerate(vcf_iterator.header.samples))
sample_index_dict = {v: k for k, v in index_sample_dict.items()}

def score_dictionary_record(record):
    supp_vec = [pos for pos, char in enumerate(record.info["SUPP_VEC"]) if char == "1"]
    sample_flag_dict = {}
    lensv = record.info["SVLEN"]
    for i in supp_vec:
            sample_name = index_sample_dict[i]
            psv = record.samples[i]["PSV"] # PSV is the support vector for this sample
            flag: int = 0
            if (psv == "100") and (lensv < 3000 or lensv > 100) : flag = 2
            elif (psv == "100") and (lensv > 3000) : flag = 1
            elif (psv == "001") : flag = 1
            elif (psv == "010") : flag = 1
            elif (psv.count("1") == 2) and (lensv < 30000): flag = 2
            elif (psv.count("1") == 3) and (lensv < 30000): flag = 3
            elif (psv.count("1") == 2) and (lensv > 30000): flag = 2
            elif (psv.count("1") == 3) and (lensv > 30000): flag = 2
            sample_flag_dict[i] = flag
    return sample_flag_dict        

for record in vcf_iterator:
    new_record = vcf_out.header.new_record(contig = record.chrom, start = record.start, stop = record.stop, alleles = record.alleles, id = record.id, qual = record.qual, filter = record.filter, info = record.info)
    key = "{}:{}:{}:{}".format(record.chrom, record.start, record.stop, record.info["SVTYPE"])
    sample_flag_dict = score_dictionary_record(record)

    for i in range(len(vcf_iterator.header.samples)):
        for format_ in record.samples[i]: new_record.samples[i][format_] = record.samples[i][format_]
        if i not in sample_flag_dict: new_record.samples[i]["SCORE"] = 0
        else: new_record.samples[i]["SCORE"] = sample_flag_dict[i]

    vcf_out.write(new_record)
