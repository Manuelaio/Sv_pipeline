import pysam
import sys

# Read the input VCF file:
delly_bcf_filepath = sys.argv[1]
delly_bcf_reader = pysam.VariantFile(delly_bcf_filepath, "r")
# New VCF file to be written:
delly_recoded_vcf_filepath = sys.argv[2]
delly_recoded_vcf_writer = pysam.VariantFile(delly_recoded_vcf_filepath, "w", header = delly_bcf_reader.header)
for record in delly_bcf_reader:
    svlen = record.stop - record.start
    if record.info["SVTYPE"] != "INS":
        record.info.__setitem__("SVLEN", svlen)
    delly_recoded_vcf_writer.write(record)
