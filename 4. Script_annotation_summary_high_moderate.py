import csv

# Input file containing all annotations from SnpSift
input_file = "annotation_summary_all.tsv"

# Output file that later will contain only HIGH or MODERATE annotations
output_file = "annotation_summary_high_moderate.tsv"

paired_columns = [
    "ANN[*].GENEID",
    "ANN[*].GENE",
    "ANN[*].IMPACT",
    "ANN[*].EFFECT",
    "ANN[*].HGVS_C",
    "ANN[*].HGVS_P"
]

with open(input_file, newline="") as infile, open(output_file, "w", newline="") as outfile:
    
    reader = csv.DictReader(infile, delimiter="\t")
    
    # Prepare the TSV writer using the same column names as the input file
    writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames, delimiter="\t")
    writer.writeheader()

    # Process each row from the input file
    for row in reader:
        
        # Store the positions of annotations with HIGH or MODERATE impact
        impacts = row["ANN[*].IMPACT"].split(",")
        keep = []
        for i, impact in enumerate(impacts):
            if impact == "HIGH" or impact == "MODERATE":
                keep.append(i)

        # Skip this row if it has no HIGH or MODERATE annotation
        if not keep:
            continue

        # Create one output row for each HIGH or MODERATE annotation
        for i in keep:
            new_row = row.copy()
            for column in paired_columns:
                values = row[column].split(",")
                new_row[column] = values[i]
            writer.writerow(new_row)

