#!/usr/bin/python 3
import argparse
import os

def read_tsv(args):
    with open(args, 'r') as file:
        lines = file.readlines()
    return lines

# to look for in maf
def parse_snvs(snvs_ls):
    mod_snvs_ls = list()
    for line in snvs_ls:
        columns = line.strip().split('\t')

        gene = columns[0] 
        if not "splice" in columns[1]:
            first_aa = columns[4].split(":")[0]
            secnd_aa = columns[8].split(":")[0]
            pos_aa = columns[1]
            mutation = f"{first_aa}{pos_aa}{secnd_aa}"
        else:
            mutation = columns[1]

        genomic_positions = columns[10].split('|')
        for genomic_position in genomic_positions:
            chromosome, position = genomic_position.split(':')[0], genomic_position.split(':')[1].split('_')[0]
            ## bed coordinate:
            # position = int(position) - 1
            mod_snvs_ls.append(f"chr{chromosome}\t{position}\t{position}\t{gene}\t{mutation}")
    
    return mod_snvs_ls

# to look for in maf
def parse_indels(indels_ls):
    mod_indels_ls = list()

    for line in indels_ls:
        columns = line.strip().split('\t')

        gene = columns[0] 
        mutation = columns[9].split(":")[0]

        genomic_positions = columns[11].split('|')
        for genomic_position in genomic_positions:
            chromosome, position = genomic_position.split(':')[0], genomic_position.split(':')[1].split('_')[0]
            mod_indels_ls.append(f"chr{chromosome}\t{position}\t{position}\t{gene}\t{mutation}")
            
    return mod_indels_ls


def maf_vs_tsv(tsv_ls, maf_ls):
    
    bed_set = set()
    # print(tsv_ls)

    for line in maf_ls:
        lines_ls = line.strip().split('\t')
        search_line = f"{lines_ls[0]}\t{lines_ls[1]}\t{lines_ls[2]}\t{lines_ls[3]}\t{lines_ls[4]}\n"
        # fewer_lines = "\t".join(lines_ls[0:5])
        # print(search_line)
        # c_change = lines_ls[5]
        if search_line in tsv_ls:
            # print(search_line)
            lines_ls[1] = str(int(lines_ls[1]) - 1)
            lines_ls = "\t".join(lines_ls)
            print(lines_ls)
            bed_set.add(lines_ls)
    bed_ls = list(bed_set)
    return bed_ls

def save_tsv(args,mod_ls):
    with open(args, 'w') as file:
        for line in mod_ls:
            file.write(f"{line}\n")

def main():
    parser = argparse.ArgumentParser(description="(1) Modify Hotspots Results V2 spreadsheet from https://www.cancerhotspots.org/#/download (run with --variant) and (2) compare maf vs tsv (run without --variant)")
    parser.add_argument("--output", required=True, help="Path to output TSV or BED file, unsorted")
    parser.add_argument("--variant",help="snv or indel")
    parser.add_argument("--tsv", required=True, help="Path to tsv generated concatenating indels and snvs or path to just indel or snv")
    parser.add_argument("--maf", help="Path to simplified maf generated with mutations found in lung cancer samples")
  
    args = parser.parse_args()

    errors = []
    if not os.path.exists(os.path.abspath(args.tsv)):
        errors.append(f"The input file {args.tsv} does not exist")
    if os.path.exists(os.path.abspath(args.output)):
        errors.append(f"The output file {args.output} already exists")
    if args.maf:
        if not os.path.exists(os.path.abspath(args.maf)):
            errors.append(f"The input file {args.maf} does not exist")
          
    if errors:
        for error in errors:
            print(error)
        exit("Exiting")
    else:
        if args.tsv and args.variant:     
            lines_ls = read_tsv(os.path.abspath(args.tsv))

            if args.variant.lower() == "snv":
                mod_lines_ls = parse_snvs(lines_ls)
            else:
                mod_lines_ls = parse_indels(lines_ls)
        elif args.tsv and args.maf and not args.variant:
            maf_ls = read_tsv(os.path.abspath(args.maf))
            tsv_ls = read_tsv(os.path.abspath(args.tsv))
            mod_lines_ls = maf_vs_tsv(tsv_ls, maf_ls)
            # print(mod_lines_ls)

        save_tsv(os.path.abspath(args.output), mod_lines_ls )

if __name__ == "__main__":
    main()
