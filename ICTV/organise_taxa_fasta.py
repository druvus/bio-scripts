import argparse
import logging
import os
import pandas as pd
from Bio import SeqIO
import re

def setup_logging():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract sequences based on ICTV Excel data and save into new FASTA files.")
    parser.add_argument("excel_file", help="Excel file with the virus data.")
    parser.add_argument("fasta_file", help="FASTA file with the sequences.")
    parser.add_argument("--output_folder", default=".", help="Folder to save the output FASTA files.")
    parser.add_argument("--segment", default="all", choices=["L", "M", "S", "all"], help="Segment to export (L, M, S, or all).")
    parser.add_argument("--family", help="Filter by family.")
    parser.add_argument("--genus", help="Filter by genus.")
    parser.add_argument("--species", help="Filter by species.")
    return parser.parse_args()

def sanitize_filename(filename):
    """Sanitize the filename by replacing spaces and special characters with underscores."""
    filename = re.sub(r'[^\w\s-]', '', filename)
    filename = re.sub(r'\s+', '_', filename)
    return filename

def extract_data_from_excel(excel_path, segment, family=None, genus=None, species=None):
    """Extract relevant data from the Excel file and filter by segment if not 'all' and other taxonomy filters."""
    df = pd.read_excel(excel_path)
    if segment != "all":
        df = df[df['Virus GENBANK accession'].str.contains(f'{segment}:')]
    if family:
        df = df[df['Family'].str.contains(family, case=False, na=False)]
    if genus:
        df = df[df['Genus'].str.contains(genus, case=False, na=False)]
    if species:
        df = df[df['Species'].str.contains(species, case=False, na=False)]
    return df[['Species', 'Virus name(s)', 'Virus GENBANK accession']].dropna()

def parse_accession_ids(accession_str, segment):
    """Parse the accession IDs from the string, handling multiple entries and ranges, filtered by segment if not 'all'."""
    acc_dict = {}
    id_pairs = re.findall(r'\b([LSM]):\s*([\w-]+(?:,[\w-]+)*)', accession_str)
    for key, group in id_pairs:
        if segment != "all" and key != segment:
            continue
        ids = []
        for part in group.split(','):
            if '-' in part:
                start, end = part.split('-')
                start_num = re.search(r'\d+$', start).group(0)
                end_num = re.search(r'\d+$', end).group(0)
                start_prefix = start[:-len(start_num)]
                for num in range(int(start_num), int(end_num) + 1):
                    ids.append(f"{start_prefix}{num}")
            else:
                ids.append(part)
        acc_dict[key] = ids
    return acc_dict

def read_fasta_sequences(fasta_path, accessions):
    """Read and store FASTA sequences that match the given accession IDs."""
    records = SeqIO.to_dict(SeqIO.parse(fasta_path, 'fasta'), key_function=lambda record: record.id.split('.')[0])
    return {acc: records.get(acc) for acc in accessions if acc in records}

def write_fasta_files(df, sequences, output_folder, segment):
    """Write sequences to new FASTA files based on DataFrame information and segment."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    for _, row in df.iterrows():
        species = sanitize_filename(row['Species'])
        virus_name = sanitize_filename(row['Virus name(s)'])
        filename_suffix = f"_{segment}.fasta" if segment != "all" else "_all.fasta"
        filename = f"{species}_{virus_name}{filename_suffix}"
        filepath = os.path.join(output_folder, filename)
        acc_dict = parse_accession_ids(row['Virus GENBANK accession'], segment)
        with open(filepath, 'w') as f:
            for key in acc_dict:
                for acc in acc_dict[key]:
                    if sequences.get(acc):
                        SeqIO.write(sequences[acc], f, 'fasta')
        logging.info(f"Written {filepath}")

def main():
    setup_logging()
    args = parse_arguments()
    df = extract_data_from_excel(args.excel_file, args.segment, args.family, args.genus, args.species)
    accessions = set(acc for row in df['Virus GENBANK accession'] for key in parse_accession_ids(row, args.segment) for acc in parse_accession_ids(row, args.segment)[key])
    sequences = read_fasta_sequences(args.fasta_file, accessions)
    write_fasta_files(df, sequences, args.output_folder, args.segment)

if __name__ == "__main__":
    main()

