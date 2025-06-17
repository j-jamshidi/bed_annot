import pandas as pd
from pyensembl import EnsemblRelease
import os
import sys

ensembl = EnsemblRelease(105)  # GRCh38/hg38

# Load MANE Select transcript IDs
def load_mane_transcripts(mane_file):
    mane_df = pd.read_csv(mane_file, sep="\t", comment="#", usecols=["Ensembl_Gene", "Ensembl_nuc"])
    mane_df.columns = ["gene", "transcript"]
    return dict(zip(mane_df["gene"], mane_df["transcript"]))

mane_transcripts = load_mane_transcripts("MANE.GRCh38.v1.4.summary.txt")

def get_gene_exon_numbers_mane(chrom, start, end):
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    genes = ensembl.genes_at_locus(contig=chrom, position=start)
    if not genes:
        return ".", ".", "."

    for gene in genes:
        gene_name = gene.gene_name
        mane_tx_id = mane_transcripts.get(gene.gene_id)

        if mane_tx_id:
            try:
                transcript = ensembl.transcript_by_id(mane_tx_id)
            except Exception as e:
                print(f"[WARN] MANE transcript {mane_tx_id} not found for gene {gene_name}: {e}")
                transcript = None
        else:
            transcript = None

        # Fallback to longest transcript if MANE not available
        if not transcript:
            if not gene.transcripts:
                continue
            transcript = max(gene.transcripts, key=lambda t: t.end - t.start)

        if not transcript or not transcript.exons:
            continue

        exon_numbers = [
            i + 1 for i, exon in enumerate(transcript.exons)
            if exon.end >= start and exon.start <= end
        ]
        if exon_numbers:
            return gene_name, f"exon{min(exon_numbers)}", f"exon{max(exon_numbers)}"

        # if transcript found but no exon matched
        return gene_name, ".", "."

    return ".", ".", "."


def annotate_bed(input_bed, output_bed):
    cols = ["chrom", "start", "end"]
    df = pd.read_csv(input_bed, sep="\t", header=None, names=cols[:3])

    annotations = []
    for _, row in df.iterrows():
        gene, exon_start, exon_end = get_gene_exon_numbers_mane(row["chrom"], row["start"], row["end"])
        annotations.append([gene, exon_start, exon_end])

    annotated_df = pd.concat([df, pd.DataFrame(annotations, columns=["gene", "first_exon", "last_exon"])], axis=1)
    annotated_df.to_csv(output_bed, sep="\t", header=False, index=False)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python geneannot.py <input.bed>")
        sys.exit(1)

    input_bed = sys.argv[1]
    bed_basename = os.path.basename(input_bed)
    bed_root, _ = os.path.splitext(bed_basename)
    output_bed = f"annot_{bed_root}.bed"

    annotate_bed(input_bed, output_bed)
    print(f"Annotation saved to {output_bed}")
