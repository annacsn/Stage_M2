import gzip
import csv
import argparse

def is_mutant(gt):
    return any(x in gt for x in ["1/1", "1|1", "1/0", "0/1", "1|0", "0|1"])

def parse_info_field(info):
    aa3to1 = {
        'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
        'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
        'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
        'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
        'Sec': 'U', 'Ter': '*'
    }
    for entry in info.split(';'):
        if entry.startswith("ANN="):
            annotations = entry[4:].split(",")
            for ann in annotations:
                fields = ann.split('|')
                if len(fields) > 9 and fields[3] == 'K13' and fields[1] == 'missense_variant':
                    protein_change = fields[10]  # ex: p.Cys580Tyr
                    if protein_change.startswith('p.'):
                        protein_change = protein_change[2:]
                    for aa3 in aa3to1:
                        if protein_change.startswith(aa3):
                            ref = aa3to1[aa3]
                            rest = protein_change[len(aa3):]
                            for aa3b in aa3to1:
                                if rest.endswith(aa3b):
                                    pos = rest[:-len(aa3b)]
                                    alt = aa3to1[aa3b]
                                    return f"{ref}{pos}{alt}"
    return None

def parse_vcf(vcf_path, sample_mutations):
    open_func = gzip.open if vcf_path.endswith('.gz') else open
    with open_func(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                samples = header[9:]
                for sample in samples:
                    if sample not in sample_mutations:
                        sample_mutations[sample] = set()
                continue

            cols = line.strip().split('\t')
            info = cols[7]
            mutation = parse_info_field(info)
            if not mutation:
                continue

            formats = cols[8].split(':')
            for i, sample in enumerate(samples):
                fields = cols[9 + i].split(':')
                fmt = dict(zip(formats, fields))
                gt = fmt.get("GT", "./.")
                if gt in ["./.", ".|."]:
                    continue
                if is_mutant(gt):
                    sample_mutations[sample].add(mutation)

def update_metadata(metadata_path, output_path, sample_mutations):
    with open(metadata_path, 'r') as meta_in, open(output_path, 'w') as meta_out:
        reader = csv.DictReader(meta_in, delimiter='\t')
        fieldnames = reader.fieldnames + ["Mutation"]
        writer = csv.DictWriter(meta_out, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in reader:
            sample = row['Sample']
            muts = sample_mutations.get(sample, set())
            row['Mutation'] = ','.join(sorted(muts)) if muts else "WT"
            writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Annoter les mutations K13 dans un fichier metadata à partir d’un VCF.")
    parser.add_argument('-v', '--vcf', required=True, help="Fichier VCF (peut être .vcf ou .vcf.gz)")
    parser.add_argument('-m', '--metadata', required=True, help="Fichier metadata d’entrée (.txt)")
    parser.add_argument('-o', '--output', required=True, help="Fichier de sortie avec les mutations ajoutées")

    args = parser.parse_args()

    sample_mutations = {}
    parse_vcf(args.vcf, sample_mutations)
    update_metadata(args.metadata, args.output, sample_mutations)

if __name__ == "__main__":
    main()