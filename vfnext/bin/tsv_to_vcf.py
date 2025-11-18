#!/usr/bin/env python3
import sys

def convert_tabular_to_vcf(tabular_path, vcf_path, sample_name):
    with open(tabular_path, "r") as inp, open(vcf_path, "w") as out:
        # VCF header
        out.write("##fileformat=VCFv4.2\n")
        out.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")
        out.write("##FILTER=<ID=FAIL,Description=\"Some of the filters failed\">\n")
        out.write("##commandline=samtools mpileup -aa -d 50000 --reference <ref> -a -B <sortedbam> | ivar variants -p <sample_id> -q <map_quality> -t 0.05\n")
        out.write("##commandline=python $projectDir/bin/tsv_to_vcf.py <sample_id>.tsv <sample_id>.ivar.vcf <sample_name>\n")
        out.write("##INFO=<ID=REF_DP,Number=1,Type=Integer,Description=\"Reference depth\">\n")
        out.write("##INFO=<ID=REF_RV,Number=1,Type=Integer,Description=\"Reference reverse strand count\">\n")
        out.write("##INFO=<ID=REF_QUAL,Number=1,Type=Float,Description=\"Reference quality\">\n")
        out.write("##INFO=<ID=ALT_DP,Number=1,Type=Integer,Description=\"Alt depth\">\n")
        out.write("##INFO=<ID=ALT_RV,Number=1,Type=Integer,Description=\"Alt reverse strand count\">\n")
        out.write("##INFO=<ID=ALT_QUAL,Number=1,Type=Float,Description=\"Alt quality\">\n")
        out.write("##INFO=<ID=ALT_FREQ,Number=1,Type=Float,Description=\"Variant allele frequency\">\n")
        out.write("##INFO=<ID=TOTAL_DP,Number=1,Type=Integer,Description=\"Total depth\">\n")
        out.write("##INFO=<ID=PVAL,Number=1,Type=Float,Description=\"Variant significance p-value\">\n")
        out.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")

        header = True
        for line in inp:
            if header:
                header = False
                continue

            fields = line.strip().split("\t")
            if len(fields) < 14:
                continue
            alt_freq = fields[10]
            # Skip variants with alt frequency < 0.5 (minor variants)
            if float(alt_freq) < 0.5:
                continue
            chrom = fields[0]
            pos = fields[1]
            ref = fields[2].upper()
            alt_raw = fields[3].upper()

            ref_dp = fields[4]
            ref_rv = fields[5]
            ref_qual = fields[6]
            alt_dp = fields[7]
            alt_rv = fields[8]
            alt_qual = fields[9]
            total_dp = fields[11]
            pval = fields[12]
            passed = fields[13]

            filter_val = "PASS" if passed == "TRUE" else "FAIL"

            # Detect and convert INDEL format
            if alt_raw.startswith("+"):  # insertion
                ins = alt_raw[1:]
                alt = ref + ins
            elif alt_raw.startswith("-"):  # deletion
                deletion = alt_raw[1:]
                alt = ref
                ref = ref + deletion
            else:  # SNV
                alt = alt_raw

            qual = "."  # no quality provided

            info = (
                f"REF_DP={ref_dp};REF_RV={ref_rv};REF_QUAL={ref_qual};"
                f"ALT_DP={alt_dp};ALT_RV={alt_rv};ALT_QUAL={alt_qual};"
                f"ALT_FREQ={alt_freq};TOTAL_DP={total_dp};PVAL={pval}"
            )

            out.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\t{filter_val}\t{info}\tGT\t1\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("\nUso:\n  python3 convert_tabular_to_vcf.py input.tsv output.vcf sample_name\n")
        sys.exit(1)

    convert_tabular_to_vcf(sys.argv[1], sys.argv[2], sys.argv[3])
