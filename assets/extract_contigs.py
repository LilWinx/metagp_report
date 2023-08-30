def read_contig_names(filename):
    with open(filename, "r") as file:
        return set(line.strip() for line in file)

def extract_match_contigs(assembly_filename, match_contigs_filename, output_filename):
    match_contigs = read_contig_names(match_contigs_filename)
    
    with open(assembly_filename, "r") as assembly_file, \
         open(output_filename, "w") as output_file:
        current_contig = None
        keep_contig = False
        
        for line in assembly_file:
            if line.startswith(">"):
                current_contig = line.strip()[1:].split()[0]
                keep_contig = current_contig in match_contigs
                if keep_contig:
                    output_file.write(line)
            elif keep_contig:
                output_file.write(line)

if __name__ == "__main__":
    assembly_filename = "assembly.fasta"  # Replace with your assembly file
    match_contigs_filename = "match_contigs.txt"  # Replace with your contig list file
    output_filename = "extracted_match_contigs.fasta"  # Replace with your desired output file
    
    extract_match_contigs(assembly_filename, match_contigs_filename, output_filename)
    print("Extraction complete.")
