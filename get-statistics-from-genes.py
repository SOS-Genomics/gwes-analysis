import argparse
import os

# Python 3.x

# OBJECTIVE
# Get statistics of genes and pair of genes from table generated by get-genes-from-positions.py

# ARGUMENTS
# --input_path       Path to read input files.
# --gene_pairs_file    Pairs of genes file:
#                       position_1, start_pos_1, end_pos_1, gene_type_1, gene_name_1,
#                       product_1, protein_id_1, note_1,
#                       position_2, start_pos_2, end_pos_2, gene_type_2, gene_name_2,
#                       product_2, protein_id_2, note_2,
#                       distance, label, score
#                       Tab-separated-value file.
# --output_path      Path to place output files.
# --output_file      Output statistics file:
#                       Tab-separated-value file. WARNING: Include inside headers for each stats table.

# OUTPUT
# 1) Output file:
# Several statistics in table format: 
# SCORE OF PAIR OF GENES HIGHLY ASSOCIATED
#   Gene 1    Gene 2  Score
# TOTAL OF PAIR OF GENES HIGHLY ASSOCIATED
#   Gene 1  Gene 2  Total
# TOTAL GENES HIGHLY ASSOCIATED (POSITION 1)
#   Gene    Total
# TOTAL GENES HIGHLY ASSOCIATED (POSITION 2)
#   Gene    Total

# RUN:
# python get-statistics-from-genes.py
# --input_path /home/user/my_path
# --gene_pairs_file genes_highly_associated.csv
# --output_path /home/user/my_path
# --output_file stats_genes_highly_associated.txt

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description='Get statistics of genes and pair of genes from table generated by get-genes-from-positions.py.')
    parser.add_argument("--input_path", dest="input_path",
                      help="Path to read input files", metavar="PATH")
    parser.add_argument("--gene_pairs_file", dest="gene_pairs_file",
                      help="Gene pairs file. Tab-separated-value file: position_1, start_pos_1, end_pos_1, gene_type_1, gene_name_1, product_1, protein_id_1, note_1, position_2, start_pos_2, end_pos_2, gene_type_2, gene_name_2, product_2, protein_id_2, note_2, distance, label, score.", metavar="FILE")
    parser.add_argument("--output_path", dest="output_path",
                      help="Path to place output files", metavar="PATH")
    parser.add_argument("--output_file", dest="output_file",
                      help="Output statistics file. Tab-separated-value file. WARNING: Include inside headers for each stats table.", metavar="FILE")

    args = parser.parse_args()

    print('-------------------------------- RUNNING --------------------------------')
    print("Get statistics of genes and pair of genes from table generated by get-genes-from-positions.py."
          "\nInput: Tab-separated: position_1, start_pos_1, end_pos_1, gene_type_1, gene_name_1, product_1, protein_id_1, note_1, position_2, start_pos_2, end_pos_2, gene_type_2, gene_name_2, product_2, protein_id_2, note_2, distance, label, score."
          "\nOutput: Several statistics in table format. Tab-separated-value file. WARNING: Include inside headers for each stats table.")

    # Print parameter values
    print('-------------------------------- PARAMETERS --------------------------------')
    print("Path to read input files: " + str(args.input_path))
    print("Gene pairs file (generated by get-genes-from-positions.py): " + str(args.gene_pairs_file))
    print("Path to place output files: " + str(args.output_path))
    print("Output statistics file: " + str(args.output_file))

    # DICTIONARIES
    hash_gene_pairs = {}
    hash_genes_pos_1 = {}
    hash_genes_pos_2 = {}
    list_gene_pairs_score = []
    print('-------------------------------- READING GENE PAIRS FILE --------------------------------')
    print('Filling out dictionaries and lists.')
    with open(os.path.join(args.input_path, args.gene_pairs_file), encoding='utf8') as iFile:
        for row in iFile:
            if row.startswith("position_1"):
                continue
            row = row.strip('\n')
            listRow = row.split('\t')
            gene_1 = listRow[4]
            gene_2 = listRow[12]
            score = listRow[18]
            if gene_1 in hash_genes_pos_1:
                hash_genes_pos_1[gene_1] += 1
            else:
                hash_genes_pos_1[gene_1] = 1

            if gene_2 in hash_genes_pos_2:
                hash_genes_pos_2[gene_2] += 1
            else:
                hash_genes_pos_2[gene_2] = 1

            gene_pair = gene_1 + "\t" + gene_2
            if gene_pair in hash_gene_pairs:
                hash_gene_pairs[gene_pair] += 1
            else:
                hash_gene_pairs[gene_pair] = 1

            list_gene_pairs_score.append((gene_pair, float(score)))
    print("File read. Dictionaries and lists done!")

    with open(os.path.join(args.output_path, args.output_file), encoding='utf8', mode="w") as oFile:
        oFile.write("-------------------------------- SCORE OF PAIR OF GENES HIGHLY ASSOCIATED --------------------------------\n")
        oFile.write("Gene 1\tGene 2\tScore\n")
        for r in sorted(list_gene_pairs_score, key=lambda x: x[1], reverse=True):
            oFile.write("{}\t{}\n".format(r[0], r[1]))
        oFile.write("-------------------------------- TOTAL OF PAIR OF GENES HIGHLY ASSOCIATED --------------------------------\n")
        oFile.write("Gene 1\tGene 2\tTotal\n")
        for k, v in sorted(hash_gene_pairs.items(), key=lambda item: item[1], reverse=True):
            oFile.write("{}\t{}\n".format(k, v))
        oFile.write("-------------------------------- TOTAL GENES HIGHLY ASSOCIATED (POSITION 1) --------------------------------\n")
        oFile.write("Gene\tTotal\n")
        for k, v in sorted(hash_genes_pos_1.items(), key=lambda item: item[1], reverse=True):
            oFile.write("{}\t{}\n".format(k, v))
        oFile.write("-------------------------------- TOTAL GENES HIGHLY ASSOCIATED (POSITION 2) --------------------------------\n")
        oFile.write("Gene\tTotal\n")
        for k, v in sorted(hash_genes_pos_2.items(), key=lambda item: item[1], reverse=True):
            oFile.write("{}\t{}\n".format(k, v))
    print('-------------------------------- DONE! --------------------------------')
