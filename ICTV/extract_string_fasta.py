import argparse
import logging
from Bio import SeqIO

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract sequences from a FASTA file based on provided search strings.")
    parser.add_argument("--fasta_file", required=True, help="The FASTA file to search through.")
    parser.add_argument("--search_strings", required=True, help="Comma-separated list of search strings to match against sequence descriptions.")
    parser.add_argument("--output_file", default="filtered_sequences.fasta", help="File to save the extracted sequences.")
    return parser.parse_args()

def filter_sequences(fasta_path, search_strings, output_file):
    search_terms = search_strings.split(',')
    logging.info(f"Searching for sequences that match any of: {search_terms}")
    sequences_found = []

    with open(fasta_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if any(search_term.lower() in record.description.lower() for search_term in search_terms):
                sequences_found.append(record)
                logging.debug(f"Matched: {record.description}")

    if sequences_found:
        SeqIO.write(sequences_found, output_file, "fasta")
        logging.info(f"Extracted {len(sequences_found)} sequences to {output_file}")
    else:
        logging.info("No matches found.")

def main():
    setup_logging()
    args = parse_arguments()
    filter_sequences(args.fasta_file, args.search_strings, args.output_file)

if __name__ == "__main__":
    main()

