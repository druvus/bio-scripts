import argparse
import regex
import logging
import gzip
import csv
from Bio import SeqIO

def count_homopolymers(seq, max_length):
    """Given a string sequence, return a dictionary counting its homopolymers"""
    # Split the string into homopolymers
    segments = filter(None, regex.findall(r'[A]+|[T]+|[C]+|[G]+', seq.upper()))
    mers_dict = {base: [0] * max_length for base in 'ATCG'}
    for segment in segments:
        ntide = segment[0]
        cnt = len(segment)
        if cnt <= max_length:
            mers_dict[ntide][cnt - 1] += 1
    return mers_dict

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Count homopolymers in a fasta file.')
    parser.add_argument('fasta', type=str, help='Fasta file (compressed or uncompressed) with sequences to count homopolymers.')
    parser.add_argument('--max_length', type=int, default=7, help='Maximum length of homopolymers to count.')
    parser.add_argument('--output', type=str, help='Output tab-delimited text file to save the results.')
    return parser.parse_args()

def open_fasta_file(filename):
    """Open fasta file, supporting both compressed (.gz) and uncompressed formats"""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def main():
    # Set up logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    args = parse_args()
    
    logging.info(f"Processing fasta file {args.fasta}")
    totals = {base: [0] * args.max_length for base in 'ATCG'}
    
    with open_fasta_file(args.fasta) as inFile:
        for record in SeqIO.parse(inFile, "fasta"):
            seq_dict = count_homopolymers(str(record.seq), args.max_length)
            logging.info(f"Processing sequence {record.name}")
            for base, counts in seq_dict.items():
                totals[base] = [sum(i) for i in zip(totals[base], counts)]
    
    results = []
    for base, counts in totals.items():
        for length in range(1, args.max_length + 1):
            nucleotide_sequence = base * length
            results.append((f"{length}{base}", nucleotide_sequence, counts[length - 1]))
    
    # Output results to tab-delimited text file if specified, otherwise print to console
    if args.output:
        logging.info(f"Writing results to {args.output}")
        with open(args.output, 'w', newline='') as txtfile:
            writer = csv.writer(txtfile, delimiter='\t')
            writer.writerow(["length", "nucleotides", "counts"])
            writer.writerows(results)
    else:
        print("length\tnucleotides\tcounts")
        for row in results:
            print("\t".join(map(str, row)))

if __name__ == '__main__':
    main()
