#PURPOSE: clean topGO enrichment files to exclude categories that have <=1 expected genes

import os
import pandas as pd

# Directory containing the enrichment CSVs
ENRICH_DIR = "/ihome/nclark/jjd108/summer_2026/rerconverge_corrected_codon/enrichment/binary_enrich/250_genes"



# The exact files to clean
FILES = [
        "RERconverge.LAURA.binary.noCOMP._GOenrichment.top250accOfAllAcc.translated.keepEarlyStops.noTraitWinsor.csv",
"RERconverge.LAURA.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv",
"RERconverge.LAURA.binary.noCOMP._GOenrichment.top250decOfAllDec.translated.keepEarlyStops.noTraitWinsor.noAsk.csv",
"RERconverge.LAURA.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv",
"RERconverge.LESS.STRICT.binary.noCOMP._GOenrichment.top250accOfAllAcc.translated.keepEarlyStops.noTraitWinsor.csv",
"RERconverge.LESS.STRICT.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv",
"RERconverge.LESS.STRICT.binary.noCOMP._GOenrichment.top250decOfAllDec.translated.keepEarlyStops.noTraitWinsor.noAsk.csv",
"RERconverge.LESS.STRICT.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv",
"RERconverge.PRIMATE.binary.noCOMP._GOenrichment.top250accOfAllAcc.translated.keepEarlyStops.noTraitWinsor.csv",
"RERconverge.PRIMATE.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv",
"RERconverge.PRIMATE.binary.noCOMP._GOenrichment.top250decOfAllDec.translated.keepEarlyStops.noTraitWinsor.noAsk.csv",
"RERconverge.PRIMATE.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv",
"RERconverge.RODENT.binary.noCOMP._GOenrichment.top250accOfAllAcc.translated.keepEarlyStops.noTraitWinsor.csv",
"RERconverge.RODENT.binary.noCOMP._GOenrichment.top250accOfAll.translated.keepEarlyStops.NoTraitWinsor.noAsk.csv",
"RERconverge.RODENT.binary.noCOMP._GOenrichment.top250decOfAllDec.translated.keepEarlyStops.noTraitWinsor.noAsk.csv",
"RERconverge.RODENT.binary.noCOMP._GOenrichment.top250decOfAll.translated.keepEarlyStops.noTraitWinsor.noAsk.csv"
        ]

for fname in FILES:
    path = os.path.join(ENRICH_DIR, fname)

    print("\n=== Processing:", fname, "===")
    df = pd.read_csv(path, sep=",")
    print("df shape:", df.shape)

    # Filter on Expected >= 1
    df["Expected"] = pd.to_numeric(df["Expected"], errors="coerce")
    df = df[df["Expected"] >= 1]
    print("new shape:", df.shape)

    # Build the output filename:
    # insert "_CLEANED.expect_greaterthan1_" right after "GOenrichment"
    out_fname = fname.replace(
        "GOenrichment.", "GOenrichment_CLEANED.expect_greaterthan1_.", 1
    )
    out_path = os.path.join(ENRICH_DIR, out_fname)

    df.to_csv(out_path, index=False)
    print("wrote:", out_path)

print("\nDone. Processed", len(FILES), "files.")
