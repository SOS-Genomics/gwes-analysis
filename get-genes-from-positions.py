import argparse
import os

# Python 3.x

# OBJECTIVE
# Get gene names from a list of pair genome positions using GenBank Feature Table file.
# Each pair includes:  position 1, position 2, distance, label, association-score

# ARGUMENTS
# --input_path       Path to read input files.
# --position_file    Pairs of positions file: position 1, position 2, distance, label, association-score. 
#                       Tab-separated-value file.
# --genebank_file    GenBank file
# --output_path      Path to place output files.
# --output_file      Output file: 
#                       position_1, start_pos_1, end_pos_1, gene_type_1, gene_name_1,
#                       product_1, protein_id_1, note_1,
#                       position_2, start_pos_2, end_pos_2, gene_type_2, gene_name_2,
#                       product_2, protein_id_2, note_2,
#                       distance, label, score
#                       Tab-separated-value file.
# --threshold       Threshold to filter out association score. Float.

# OUTPUT
# 1) Output file:
# position_1, start_pos_1, end_pos_1, gene_type_1, gene_name_1, product_1, protein_id_1, note_1,
# position_2, start_pos_2, end_pos_2, gene_type_2, gene_name_2, product_2, protein_id_2, note_2,
# distance, label, score
# 2) Missing positions:
# filename = output_file + _missing_pos.txt

# RUN:
# python get-genes-from-positions.py
# --input_path /home/user/my_path
# --position_file positions_highly_associated.txt
# --genebank_file sequence-complete-genome.gb
# --output_path /home/user/my_path
# --output_file genes_highly_associated.csv
# --threshold 1.00

if __name__ == "__main__":
    # Arguments
    parser = argparse.ArgumentParser(description='Get gene names from a list of pair genome positions using GenBank Feature Table file.')
    parser.add_argument("--input_path", dest="input_path", required=True,
                      help="Path to read input files", metavar="PATH")
    parser.add_argument("--position_file", dest="position_file", required=True,
                      help="Tab-separated: position 1, position 2, distance, label, association-score", metavar="FILE")
    parser.add_argument("--genebank_file", dest="genebank_file", required=True,
                      help="GenBank file", metavar="FILE")
    parser.add_argument("--output_path", dest="output_path", required=True,
                      help="Path to place output files", metavar="PATH")
    parser.add_argument("--output_file", dest="output_file", required=True,
                      help="Tab-separated: position_1, start_pos_1, end_pos_1, gene_type_1, gene_name_1, product_1, protein_id_1, note_1, position_2, start_pos_2, end_pos_2, gene_type_2, gene_name_2, product_2, protein_id_2, note_2, distance, label, score.", metavar="FILE")
    parser.add_argument("--threshold", dest="threshold", default=0, type=float,
                      help="Threshold to filter out association score (float).", metavar="NUMERIC")

    args = parser.parse_args()

    print('-------------------------------- RUNNING --------------------------------')
    print("Get gene names from a list of pair genome positions using GenBank Feature Table file."
          "\nInput: Tab-separated: position 1, position 2, distance, label, association-score."
          "\nInput: GenBank file."
          "\nOutput: Tab-separated: position_1, start_pos_1, end_pos_1, gene_type_1, gene_name_1, product_1, protein_id_1, note_1, position_2, start_pos_2, end_pos_2, gene_type_2, gene_name_2, product_2, protein_id_2, note_2, distance, label, score."
          "\nMissing positions: output_file + _missing_pos.txt")

    # Print parameter values
    print('-------------------------------- PARAMETERS --------------------------------')
    print("Path to read input files: " + str(args.input_path))
    print("File to read pairs of positions: " + str(args.position_file))
    print("GenBank file: " + str(args.genebank_file))
    print("Path to place output files: " + str(args.output_path))
    print("File to place genes and descriptions of pairs of positions: " + str(args.output_file))
    print("Threshold to filter out association score: " + str(args.threshold))

    # DICTIONARY: gene names and gene descriptions
    hash_genes_report = {}

    print('-------------------------------- READING POSITION FILE --------------------------------')
    # Read file with pair of positions
    input_positions = []
    with open(os.path.join(args.input_path, args.position_file), encoding="utf8") as pFile:
        for row in pFile:
            row = row.strip('\n')
            if row != "":
                listRow = row.split("\t")
                association_score = float(listRow[4])
                if association_score >= args.threshold:
                    input_positions.append(row)
    print("File read! {} pairs of positions.".format(len(input_positions)))

    # Create dictionary of unique positions and list of position tuples
    positions = []
    for row in input_positions:
        listRow = row.split('\t')
        positions.append(listRow)
        pairPositions = []
        pairPositions.append(listRow[0])
        pairPositions.append(listRow[1])
        for position in pairPositions:
            int_position = int(position)
            if position not in hash_genes_report:
                hash_genes_report[position] = {"gene_type": "",
                                          "gene_name": "",
                                          "start_pos": "",
                                          "end_pos": "",
                                          "product": "",
                                          "protein_id": "",
                                          "note": ""}
    # print("Length position pairs (list positions): {}".format(len(positions)))
    # print("Length unique positions (hash_genes_report): {}".format(len(hash_genes_report)))

    print('-------------------------------- LOOKING FOR GENE AND GENE DESCRIPTION IN GENEBANK FILE --------------------------------')

    # LIST: missing_positions
    missing_positions = []
    # READ GenBank file
    with open(os.path.join(args.input_path, args.genebank_file), encoding='utf8') as iFile:
        rows = iFile.readlines()
    print("GeneBank file read!")
    # Row-by-row look-up
    for position in hash_genes_report.keys():
        # print("Searching position {}...".format(position))
        i = 1
        while i < len(rows):
            # print("i: {}".format(i))
            row = rows[i].rstrip('\n')
            if row == "":
                i += 1
                continue
            # print("row: {}".format(row))
            listRow = row.split('\t')
            # Check for empty position 1
            if listRow[0] == "":
                start_pos = 0
            else:
                start_pos = int(listRow[0].strip("<"))
            # print("start_pos: {}".format(start_pos))

            # Check for empty position 2
            if listRow[1] == "":
                end_pos = 0
            else:
                end_pos = int(listRow[1].strip())
            # print("end_pos: {}".format(end_pos))

            # Check for reversed start and end position in GeneBank file
            if start_pos > end_pos:
                x = start_pos
                start_pos = end_pos
                end_pos = x

            # Check for position between start and end position
            if (start_pos <= int(position) and end_pos >= int(position)):
                # Get type: gene o CDS
                gtype = listRow[2]
                # print("Type: {}".format(gtype))
                # print("Position {} found! Start: {} End: {} Type: {}".format(position, start_pos, end_pos, gtype))
                if gtype == "gene":
                    # Next row for getting gene name
                    i += 1
                    row = rows[i].rstrip('\n')
                    listRow = row.split('\t')
                    gene_type = listRow[3]
                    gene_name = listRow[4]
                    hash_genes_report[position]["gene_type"] = gene_type
                    hash_genes_report[position]["gene_name"] = gene_name
                    hash_genes_report[position]["start_pos"] = start_pos
                    hash_genes_report[position]["end_pos"] = end_pos
                elif gtype == "CDS":
                    # Next rows for getting product, protein_id and note
                    # PRODUCT
                    i += 1
                    row = rows[i].rstrip('\n')
                    listRow = row.split('\t')
                    product = listRow[4]
                    hash_genes_report[position]["product"] = product
                    # PROTEIN_ID
                    i += 2
                    row = rows[i].rstrip('\n')
                    listRow = row.split('\t')
                    protein_id = listRow[4]
                    hash_genes_report[position]["protein_id"] = protein_id
                    # NOTE
                    i += 1
                    row = rows[i].rstrip('\n')
                    listRow = row.split('\t')
                    note = listRow[4]
                    hash_genes_report[position]["note"] = note
                    break
                else:
                    missing_positions.append(gtype + "\t" + position)
            i += 1
        if i >= len(rows):
            if position not in missing_positions:
                missing_positions.append("Missing\t" + position)

    print('-------------------------------- SAVING OUTPUT DATA --------------------------------')
    with open(os.path.join(args.output_path, args.output_file), encoding='utf8', mode="w") as oFile:
        oFile.write("position_1\tstart_pos_1\tend_pos_1\tgene_type_1\tgene_name_1\tproduct_1\tprotein_id_1\tnote_1\t")
        oFile.write("position_2\tstart_pos_2\tend_pos_2\tgene_type_2\tgene_name_2\tproduct_2\tprotein_id_2\tnote_2\t")
        oFile.write("distance\tlabel\tscore\n")
        # print(positions)
        for pairs in positions:
            # print(pairs)
            wrow = ""
            pos_1 = pairs[0]
            v = hash_genes_report[pos_1]
            wrow += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(pos_1, v["start_pos"], v["end_pos"], v["gene_type"],
                                                              v["gene_name"], v["product"], v["protein_id"], v["note"])
            pos_2 = pairs[1]
            v = hash_genes_report[pos_2]
            wrow += "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(pos_2, v["start_pos"], v["end_pos"], v["gene_type"], v["gene_name"], v["product"], v["protein_id"], v["note"])
            distance = pairs[2]
            label = pairs[3]
            score = pairs[4]
            wrow += "{}\t{}\t{}\n".format(distance, label, score)
            oFile.write("{}".format(wrow))
    print("Output saved in {}".format(os.path.join(args.output_path, args.output_file)))

    # PRINT MISSING POSITIONS
    print('-------------------------------- SAVE MISSING POSITIONS --------------------------------')
    print("Missing positions: {}".format(len(missing_positions)))
    print("Saving missing positions in: {}".format(os.path.join(args.output_path, args.output_file + "_missing_pos.txt")))
    with open(os.path.join(args.output_path, args.output_file + "_missing_pos.txt"), mode="w", encoding="utf8") as oFile:
        oFile.write('-------------------------------- MISSING POSITIONS --------------------------------\n')
        oFile.write("Missing positions: {}\n".format(len(missing_positions)))
        for mp in sorted(missing_positions):
            oFile.write("{}\n".format(mp))
    print('-------------------------------- DONE! --------------------------------')
