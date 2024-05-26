import argparse
import logging
from collections import defaultdict

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Parse and split accession data.")
    parser.add_argument("input_file", help="The input file with accession data.")
    parser.add_argument("--output_prefix", default="segment", help="Output prefix for the files.")
    return parser.parse_args()

def parse_accessions(file_path):
    accessions = defaultdict(list)
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(';')
            for part in parts:
                if ':' in part:
                    key, value = part.strip().split(':')
                    accessions[key.strip()].append(value.strip())
    return accessions

def write_outputs(accessions, output_prefix):
    for key in accessions:
        with open(f"{output_prefix}_{key}.txt", 'w') as file:
            file.writelines(f"{item}\n" for item in accessions[key])
        logging.info(f"Output for {key} written to {output_prefix}_{key}.txt")

def main():
    setup_logging()
    args = parse_arguments()
    logging.info(f"Reading data from {args.input_file}")
    accessions = parse_accessions(args.input_file)
    write_outputs(accessions, args.output_prefix)

if __name__ == "__main__":
    main()

