#!/usr/bin/env python3
"""Build K=3/K=4 dominant-Q population lists for pi/FST analysis.

The PCA figures in the GWAS directory were colored from ADMIXTURE Q files.
For reproducible pi/FST groups, each sample is assigned to the population
with the largest Q value for K=3 and K=4.
"""

import itertools
from pathlib import Path
from typing import Dict, List, Tuple


BASE = Path("${ANALYSIS_DIR}")
GWAS = BASE / "05_GWAS/00_analysis/02_genomeDB"
STRUCTURE_DIR = GWAS / "08_structure"
PCA_DIR = GWAS / "07_PCA"
OLD_GROUP_DIR = Path("${DATA_DIR2}/projects/1-oil_palm/06-population/2-pi-Fst")
OUTDIR = BASE / "22_answer_reviews/00_ms/03_V3/04_figure4/K3_K4_pi_fst_pca_groups"

FAM_FILE = STRUCTURE_DIR / "all.fam"
PCA_FILE = PCA_DIR / "PCA_out.eigenvec"
EIGENVAL_FILE = PCA_DIR / "PCA_out.eigenval"

OLD_GROUPS = ["AFR", "HHG", "IDB", "SA-EG", "SEA-A", "SEA-B"]
K_VALUES = [3, 4]
Q_THRESHOLD = 0.7


def read_fam(path):
    samples = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError(f"Malformed fam line in {path}: {line!r}")
            samples.append(fields[1])
    return samples


def read_q(path, k, samples):
    rows = []
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            values = [float(x) for x in line.split()]
            if len(values) != k:
                raise ValueError(f"{path} has {len(values)} columns, expected K={k}")
            rows.append(values)
    if len(rows) != len(samples):
        raise ValueError(f"{path} has {len(rows)} rows, fam has {len(samples)} samples")
    return rows


def read_pca(path):
    pca = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 5:
                raise ValueError(f"Malformed PCA line in {path}: {line!r}")
            pca[fields[1]] = fields[2:]
    return pca


def read_variance(path):
    vals = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line:
                vals.append(float(line))
    total = sum(vals)
    return [100.0 * x / total for x in vals] if total else []


def read_old_group_map(group_dir):
    mapping = {}
    for group in OLD_GROUPS:
        path = group_dir / f"{group}.txt"
        if not path.exists():
            continue
        with path.open() as handle:
            for line in handle:
                sample = line.strip()
                if sample:
                    mapping[sample] = group
    return mapping


def write_table(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("\t".join(header) + "\n")
        for row in rows:
            handle.write("\t".join(str(x) for x in row) + "\n")


def format_float(value):
    return f"{value:.6f}"


def build_for_k(
    k,
    samples,
    pca,
    old_group,
):
    q_file = STRUCTURE_DIR / f"all.{k}.Q"
    q_rows = read_q(q_file, k, samples)

    k_label = f"K{k}"
    group_dir = OUTDIR / "groups" / k_label
    group_dir.mkdir(parents=True, exist_ok=True)

    assignments = []
    summary_rows = []
    sample_lists = {f"{k_label}_Pop{i}": [] for i in range(1, k + 1)}

    for sample, q_values in zip(samples, q_rows):
        max_idx = max(range(k), key=lambda i: q_values[i])
        max_q = q_values[max_idx]
        second_q = sorted(q_values, reverse=True)[1] if k > 1 else 0.0
        pop = f"{k_label}_Pop{max_idx + 1}"
        threshold_group = pop if max_q >= Q_THRESHOLD else "Admixed"
        sample_lists[pop].append(sample)
        pc_values = pca.get(sample, ["NA", "NA", "NA"])
        row = [
            sample,
            k_label,
            pop,
            threshold_group,
            format_float(max_q),
            format_float(second_q),
            old_group.get(sample, "NA"),
            pc_values[0] if len(pc_values) > 0 else "NA",
            pc_values[1] if len(pc_values) > 1 else "NA",
            pc_values[2] if len(pc_values) > 2 else "NA",
        ]
        row.extend(format_float(x) for x in q_values)
        assignments.append(row)

    for pop, pop_samples in sample_lists.items():
        with (group_dir / f"{pop}.txt").open("w") as handle:
            handle.write("\n".join(pop_samples) + "\n")

        max_q_values = [
            float(row[4])
            for row in assignments
            if row[2] == pop
        ]
        pure_count = sum(1 for row in assignments if row[2] == pop and row[3] != "Admixed")
        summary_rows.append(
            [
                k_label,
                pop,
                len(pop_samples),
                pure_count,
                len(pop_samples) - pure_count,
                format_float(min(max_q_values)),
                format_float(sum(max_q_values) / len(max_q_values)),
                format_float(max(max_q_values)),
            ]
        )

    assignment_header = [
        "sample",
        "K",
        "dominant_group",
        f"threshold_{Q_THRESHOLD:g}_group",
        "max_Q",
        "second_Q",
        "old_6group",
        "PC1",
        "PC2",
        "PC3",
    ] + [f"Q{i}" for i in range(1, k + 1)]
    write_table(OUTDIR / "metadata" / f"{k_label}_dominant_assignments.tsv", assignment_header, assignments)

    return list(sample_lists), assignments, summary_rows


def build_tasks(k_to_groups):
    rows = []
    for k_label in sorted(k_to_groups):
        groups = k_to_groups[k_label]
        for group in groups:
            group_file = OUTDIR / "groups" / k_label / f"{group}.txt"
            out_prefix = OUTDIR / "vcftools" / k_label / f"{group}_100kb.pi"
            rows.append(["PI", k_label, group, "NA", group_file, "NA", out_prefix])
        for pop1, pop2 in itertools.combinations(groups, 2):
            pop1_file = OUTDIR / "groups" / k_label / f"{pop1}.txt"
            pop2_file = OUTDIR / "groups" / k_label / f"{pop2}.txt"
            out_prefix = OUTDIR / "vcftools" / k_label / f"{pop1}_{pop2}_100kb_fst"
            rows.append(["FST", k_label, pop1, pop2, pop1_file, pop2_file, out_prefix])
    write_table(
        OUTDIR / "metadata" / "vcftools_tasks.tsv",
        ["task", "K", "pop1", "pop2", "pop1_file", "pop2_file", "out_prefix"],
        rows,
    )


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "metadata").mkdir(exist_ok=True)
    (OUTDIR / "vcftools").mkdir(exist_ok=True)
    (OUTDIR / "figures").mkdir(exist_ok=True)
    (OUTDIR / "logs").mkdir(exist_ok=True)

    samples = read_fam(FAM_FILE)
    pca = read_pca(PCA_FILE)
    old_group = read_old_group_map(OLD_GROUP_DIR)
    variance = read_variance(EIGENVAL_FILE)

    k_to_groups = {}
    all_summary = []
    all_assignments = []

    for k in K_VALUES:
        groups, assignments, summary = build_for_k(k, samples, pca, old_group)
        k_label = f"K{k}"
        k_to_groups[k_label] = groups
        all_summary.extend(summary)
        for row in assignments:
            all_assignments.append(row[:10])

    write_table(
        OUTDIR / "metadata" / "K3_K4_group_summary.tsv",
        ["K", "group", "n_samples", f"n_maxQ_ge_{Q_THRESHOLD:g}", f"n_maxQ_lt_{Q_THRESHOLD:g}", "min_maxQ", "mean_maxQ", "max_maxQ"],
        all_summary,
    )
    write_table(
        OUTDIR / "metadata" / "K3_K4_assignments_minimal.tsv",
        ["sample", "K", "dominant_group", f"threshold_{Q_THRESHOLD:g}_group", "max_Q", "second_Q", "old_6group", "PC1", "PC2", "PC3"],
        all_assignments,
    )
    build_tasks(k_to_groups)

    with (OUTDIR / "metadata" / "input_sources.txt").open("w") as handle:
        handle.write(f"FAM\t{FAM_FILE}\n")
        handle.write(f"PCA\t{PCA_FILE}\n")
        handle.write(f"EIGENVAL\t{EIGENVAL_FILE}\n")
        handle.write(f"STRUCTURE_DIR\t{STRUCTURE_DIR}\n")
        handle.write(f"GROUP_RULE\tdominant ADMIXTURE Q, matching PCA colored-by-structure plots\n")
        handle.write("STRICT_ASSIGNMENT\tEach sample is assigned to exactly one Pop within each K by argmax(Q); no samples are dropped for pi/FST.\n")
        handle.write(f"Q_THRESHOLD_AUDIT\t{Q_THRESHOLD:g}; threshold groups are reported but not used for vcftools\n")
        if variance:
            handle.write("PCA_VARIANCE_PERCENT\t" + "\t".join(format_float(x) for x in variance[:10]) + "\n")

    print(f"Wrote K=3/K=4 grouping files under {OUTDIR}")
    print(f"Samples: {len(samples)}")
    for row in all_summary:
        print("\t".join(str(x) for x in row[:4]))


if __name__ == "__main__":
    main()
