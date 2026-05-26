"""
Bioinformatics File Format Parser
===================================
Demonstrates reading and summarizing the core file formats used in genomics.

Formats covered:
  .fasta  — nucleotide / protein sequences
  .fastq  — sequencing reads with Phred quality scores
  .vcf    — variant calls (SNPs, indels)
  .gff3   — genome feature annotations
  .bed    — genomic interval regions
  .pdb    — protein 3D structure records
  .nwk    — Newick phylogenetic trees

Run from the repo root:
    .venv/Scripts/python utils/format_parser.py
"""

import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


# ── FASTA ──────────────────────────────────────────────────────────────────────

def parse_fasta(filepath):
    records = {}
    current_id = None
    seq_parts = []
    with open(filepath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id:
                    records[current_id] = "".join(seq_parts)
                current_id = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if current_id:
        records[current_id] = "".join(seq_parts)
    return records


def summarize_fasta(filepath):
    records = parse_fasta(filepath)
    print(f"\n  FASTA  ->  {filepath.name}  ({len(records)} sequences)")
    for seq_id, seq in records.items():
        length = len(seq)
        gc = (seq.upper().count("G") + seq.upper().count("C")) / length * 100 if length else 0
        print(f"    {seq_id:<40s}  len={length:>6}  GC={gc:5.1f}%")


# ── FASTQ ──────────────────────────────────────────────────────────────────────

def parse_fastq(filepath):
    reads = []
    with open(filepath) as fh:
        while True:
            header = fh.readline().strip()
            if not header:
                break
            seq = fh.readline().strip()
            fh.readline()  # '+'
            qual = fh.readline().strip()
            reads.append((header[1:], seq, qual))
    return reads


def phred_scores(qual_string):
    return [ord(c) - 33 for c in qual_string]


def summarize_fastq(filepath):
    reads = parse_fastq(filepath)
    print(f"\n  FASTQ  ->  {filepath.name}  ({len(reads)} reads)")
    for name, seq, qual in reads:
        scores = phred_scores(qual)
        mean_q = sum(scores) / len(scores)
        print(f"    {name:<35s}  len={len(seq):>4}  mean_Q={mean_q:5.1f}")


# ── VCF ────────────────────────────────────────────────────────────────────────

def parse_vcf(filepath):
    variants = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = cols[:8]
            info_dict = dict(
                field.split("=", 1) if "=" in field else (field, True)
                for field in info.split(";")
            )
            variants.append({
                "chrom": chrom, "pos": int(pos), "id": vid,
                "ref": ref, "alt": alt, "qual": qual, "filter": filt,
                "info": info_dict,
            })
    return variants


def summarize_vcf(filepath):
    variants = parse_vcf(filepath)
    print(f"\n  VCF    ->  {filepath.name}  ({len(variants)} variants)")
    for v in variants:
        gene = v["info"].get("GENE", "?")
        csq = v["info"].get("CONSEQUENCE", "?")
        clnsig = v["info"].get("CLNSIG", "?")
        print(f"    {v['chrom']:6s}:{v['pos']:<12d}  {v['id']:<15s}  {v['ref']:>4s}->{v['alt']:<6s}  {gene:<8s}  {csq:<30s}  {clnsig}")


# ── GFF3 ───────────────────────────────────────────────────────────────────────

def parse_gff3(filepath):
    features = []
    with open(filepath) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = cols
            attr_dict = {}
            for field in attrs.split(";"):
                if "=" in field:
                    k, v = field.split("=", 1)
                    attr_dict[k.strip()] = v.strip()
            features.append({
                "seqid": seqid, "type": ftype,
                "start": int(start), "end": int(end),
                "strand": strand, "attrs": attr_dict,
            })
    return features


def summarize_gff3(filepath):
    features = parse_gff3(filepath)
    from collections import Counter
    type_counts = Counter(f["type"] for f in features)
    print(f"\n  GFF3   ->  {filepath.name}  ({len(features)} features)")
    for ftype, count in sorted(type_counts.items()):
        print(f"    {ftype:<20s}  {count:>3} features")
    genes = [f for f in features if f["type"] == "gene"]
    for g in genes:
        name = g["attrs"].get("Name", g["attrs"].get("ID", "?"))
        span = g["end"] - g["start"]
        print(f"    gene: {name:<12s}  {g['seqid']}:{g['start']}-{g['end']}  ({span:,} bp)  strand={g['strand']}")


# ── BED ────────────────────────────────────────────────────────────────────────

def parse_bed(filepath):
    regions = []
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            chrom, start, end = cols[0], int(cols[1]), int(cols[2])
            name = cols[3] if len(cols) > 3 else "."
            regions.append({"chrom": chrom, "start": start, "end": end, "name": name})
    return regions


def summarize_bed(filepath):
    regions = parse_bed(filepath)
    total_bp = sum(r["end"] - r["start"] for r in regions)
    print(f"\n  BED    ->  {filepath.name}  ({len(regions)} regions, {total_bp:,} bp total)")
    for r in regions:
        size = r["end"] - r["start"]
        print(f"    {r['chrom']:6s}:{r['start']:<12d}-{r['end']:<12d}  {size:>7,} bp  {r['name']}")


# ── PDB ────────────────────────────────────────────────────────────────────────

def parse_pdb(filepath):
    atoms = []
    ssbonds = []
    seqres = {}
    with open(filepath) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec in ("ATOM", "HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21].strip()
                res_seq = line[22:26].strip()
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                atoms.append({"type": rec, "atom": atom_name, "res": res_name,
                               "chain": chain, "res_seq": res_seq, "xyz": (x, y, z)})
            elif rec == "SSBOND":
                ssbonds.append(line.strip())
            elif rec == "SEQRES":
                chain = line[11].strip()
                residues = line[19:].split()
                seqres.setdefault(chain, []).extend(residues)
    return atoms, ssbonds, seqres


def summarize_pdb(filepath):
    atoms, ssbonds, seqres = parse_pdb(filepath)
    chains = sorted(set(a["chain"] for a in atoms if a["chain"]))
    print(f"\n  PDB    ->  {filepath.name}")
    print(f"    Chains: {chains}  |  Total atoms: {len(atoms)}  |  Disulfide bonds: {len(ssbonds)}")
    for chain, residues in seqres.items():
        print(f"    Chain {chain}: {len(residues)} residues — {' '.join(residues[:8])} ...")


# ── NEWICK ─────────────────────────────────────────────────────────────────────

def count_newick_taxa(newick_str):
    return re.findall(r'[A-Za-z][A-Za-z0-9_]+(?=[:,()\s])', newick_str)


def summarize_newick(filepath):
    with open(filepath) as fh:
        newick = fh.read().strip()
    taxa = count_newick_taxa(newick)
    print(f"\n  NEWICK ->  {filepath.name}  ({len(taxa)} taxa)")
    for t in taxa:
        print(f"    {t}")
    # extract branch lengths
    lengths = [float(x) for x in re.findall(r':(\d+\.\d+)', newick)]
    if lengths:
        print(f"    Branch lengths: min={min(lengths):.4f}  max={max(lengths):.4f}  mean={sum(lengths)/len(lengths):.4f}")


# ── MAIN ───────────────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("  Bioinformatics File Format Demo")
    print("=" * 70)

    summarize_fasta(ROOT / "04_protein_analysis/data/proteins.fasta")
    summarize_fasta(ROOT / "04_protein_analysis/data/gene_sequences.fasta")
    summarize_fastq(ROOT / "01_rna_seq/data/sample_reads.fastq")
    summarize_vcf(ROOT / "03_population_genetics/data/variants.vcf")
    summarize_gff3(ROOT / "01_rna_seq/data/annotation.gff3")
    summarize_bed(ROOT / "03_population_genetics/data/regions.bed")
    summarize_pdb(ROOT / "04_protein_analysis/data/insulin.pdb")
    summarize_newick(ROOT / "04_protein_analysis/data/brca2_phylogeny.nwk")

    print("\n" + "=" * 70)
    print("  All formats parsed successfully.")
    print("=" * 70)


if __name__ == "__main__":
    main()
