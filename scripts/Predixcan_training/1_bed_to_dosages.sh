#!/usr/bin/env bash

BUILD="b37"
BIM="predixcan_training_data.bim"
TRAW="predixcan_training_data_dosage_trans.traw"
OUT_DIR="processed/"

mkdir -p "$OUT_DIR"

awk -v build="$BUILD" -v out_dir="$OUT_DIR" '
    BEGIN { FS="\t"; OFS="\t" }

    NR == FNR {
        chr=$1; snp_id=$2; pos=$4; alt=$5; ref=$6
        varid = chr "_" pos "_" ref "_" alt "_" build
        annot_file = out_dir "/snp_annot.chr" chr ".txt"
        if (!seen_chr[chr]++) {
            print "chromosome", "pos", "varID", "ref_vcf", "alt_vcf", "rsid" > annot_file
        }
        print chr, pos, varid, ref, alt, snp_id > annot_file
        map_varid[snp_id] = varid
        next
    }

    FNR == 1 {
        header = "varID"
        for (i=7; i<=NF; i++) {
            # Take only the part before the first underscore
            split($i, a, "_")
            header = header "\t" a[1]
        }
        next
    }

    {
        chr = $1
        rsid = $2
        varid = map_varid[rsid]
        if (varid == "") varid = rsid

        geno_file = out_dir "/genotype.chr" chr ".txt"
        if (!seen_geno_chr[chr]++) {
            print header > geno_file
        }

        row = varid
        for (i=7; i<=NF; i++) {
            row = row "\t" $i
        }
        print row > geno_file
    }
' "$BIM" "$TRAW"