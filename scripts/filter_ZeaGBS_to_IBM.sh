#!/bin/bash

XLSX="./data/AllZeaGBSv2.7_publicSamples_metadata20140411.xlsx"
FULL_GENO_FILE="./data/ZeaGBSv27_publicSamples_raw_AGPv4-181023.vcf.gz"  # set this

# Index if there isn't one
if [ ! -f "${FULL_GENO_FILE}.tbi" ] && [ ! -f "${FULL_GENO_FILE}.csi" ]; then
    echo "No index found for $FULL_GENO_FILE — building tabix index..."
    tabix -p vcf "$FULL_GENO_FILE"
fi

# Convert xlsx to tsv, filter rows, build sample string
SAMPLES_STRING=$(python3 - <<'EOF'
import openpyxl, sys

wb = openpyxl.load_workbook(sys.argv[1] if len(sys.argv) > 1 else "./data/AllZeaGBSv2.7_publicSamples_metadata20140411.xlsx", read_only=True, data_only=True)
ws = wb.active

headers = [cell.value for cell in next(ws.iter_rows(min_row=1, max_row=1))]
germ_idx = headers.index("GermplasmSet")
dna_idx  = headers.index("DNASample")
lib_idx  = headers.index("LibraryPrepID")

samples = []
for row in ws.iter_rows(min_row=2, values_only=True):
    germplasm = row[germ_idx]
    dna       = row[dna_idx]
    lib       = row[lib_idx]
    if germplasm == "IBM" or dna == "Mo17":
        samples.append(f"{dna}:{lib}")

print(",".join(samples))
EOF
)

echo "Sample string: $SAMPLES_STRING"

bcftools view -s "$SAMPLES_STRING" -Oz -o test_ZeaGBSv27_IBM_raw_AGPv4.vcf.gz "$FULL_GENO_FILE"