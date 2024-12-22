import argparse
import sys
import pandas as pd
from Bio import SeqIO
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Delete gap-causing sequences from a FASTA file based on an Excel report.")
    parser.add_argument("-r", "--report", required=True, help="Path to the Excel gap-causers report.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to the original FASTA file.")
    parser.add_argument("-k", "--keep_fasta", required=False, default="filtered_sequences.fasta",
                        help="Path to save the FASTA file with remaining sequences (default: filtered_sequences.fasta).")
    parser.add_argument("-d", "--deleted_fasta", required=False, default="deleted_sequences.fasta",
                        help="Path to save the FASTA file with deleted sequences (default: deleted_sequences.fasta).")
    return parser.parse_args()

def extract_gap_causers(report_path):
    """
    Extracts a set of truncated sequence IDs that are identified as gap-causers from the Excel report.
    
    Args:
        report_path (str): Path to the Excel report.
    
    Returns:
        set: A set of truncated sequence IDs causing gaps.
    """
    try:
        excel_file = pd.ExcelFile(report_path)
    except Exception as e:
        print(f"Error reading Excel report: {e}")
        sys.exit(1)
    
    gap_causers = set()

    for sheet_name in excel_file.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        
        if sheet_name.startswith('Small_Gap_Causers'):
            # Extract all Gap_Causers_Sequence_X columns
            seq_columns = [col for col in df.columns if col.startswith('Gap_Causers_Sequence_')]
            for col in seq_columns:
                sequences = df[col].dropna().astype(str).str.strip()
                gap_causers.update(sequences.tolist())
        
        elif sheet_name.startswith('Large_Gap_Causers'):
            # Extract Gap_Causers_Sequences_Combined column
            if 'Gap_Causers_Sequences_Combined' in df.columns:
                combined_sequences = df['Gap_Causers_Sequences_Combined'].dropna().astype(str)
                for entry in combined_sequences:
                    # Split by comma and strip whitespaces
                    sequences = [seq.strip() for seq in entry.split(',') if seq.strip()]
                    gap_causers.update(sequences)
    
    return gap_causers

def get_user_selected_positions(gap_positions):
    """
    Prompts the user to select gap positions to delete sequences from.
    
    Args:
        gap_positions (list): List of available gap positions.
    
    Returns:
        set: Set of selected gap positions.
    """
    print("\nAvailable Gap Positions:")
    print("------------------------")
    for pos in gap_positions:
        print(f"Position: {pos}")
    print("------------------------")
    print(f"Total gap positions: {len(gap_positions)}\n")
    
    print("Select gap positions to delete sequences from.")
    print("You can enter individual positions separated by commas (e.g., 100,200,300) ")
    print("or specify ranges using hyphens (e.g., 100-150,200,300-350).")
    
    user_input = input("Enter gap positions to delete: ").strip()
    
    if not user_input:
        print("No input provided. Exiting.")
        sys.exit(1)
    
    selected_positions = set()
    parts = user_input.split(',')
    for part in parts:
        part = part.strip()
        if '-' in part:
            try:
                start, end = map(int, part.split('-'))
                if start > end:
                    start, end = end, start
                selected_positions.update(range(start, end + 1))
            except ValueError:
                print(f"Invalid range input: '{part}'. Skipping.")
        else:
            try:
                pos = int(part)
                selected_positions.add(pos)
            except ValueError:
                print(f"Invalid position input: '{part}'. Skipping.")
    
    # Validate selected positions
    valid_positions = selected_positions & set(gap_positions)
    invalid_positions = selected_positions - set(gap_positions)
    
    if invalid_positions:
        print("\nWarning: The following positions are not present in the report and will be ignored:")
        print(", ".join(map(str, sorted(invalid_positions))))
    
    if not valid_positions:
        print("No valid positions selected. Exiting.")
        sys.exit(1)
    
    print(f"\nSelected {len(valid_positions)} position(s) for deletion.\n")
    return valid_positions

def collect_sequences_to_delete(selected_positions, report_path):
    """
    Collects all truncated sequence IDs associated with the selected gap positions.
    
    Args:
        selected_positions (set): Set of selected gap positions.
        report_path (str): Path to the Excel gap-causers report.
    
    Returns:
        set: Set of truncated sequence IDs to delete.
    """
    try:
        excel_file = pd.ExcelFile(report_path)
    except Exception as e:
        print(f"Error reading Excel report: {e}")
        sys.exit(1)
    
    sequences_to_delete = set()

    for sheet_name in excel_file.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        
        if sheet_name.startswith('Small_Gap_Causers'):
            # Extract all Gap_Causers_Sequence_X columns
            seq_columns = [col for col in df.columns if col.startswith('Gap_Causers_Sequence_')]
            for _, row in df.iterrows():
                position = row['Position']
                if position in selected_positions:
                    for col in seq_columns:
                        seq_id = row[col]
                        if pd.notna(seq_id):
                            sequences_to_delete.add(str(seq_id).strip())
        
        elif sheet_name.startswith('Large_Gap_Causers'):
            # Extract Gap_Causers_Sequences_Combined column
            if 'Gap_Causers_Sequences_Combined' in df.columns:
                for _, row in df.iterrows():
                    position = row['Position']
                    if position in selected_positions:
                        combined_seqs = row['Gap_Causers_Sequences_Combined']
                        if pd.notna(combined_seqs):
                            sequences = [seq.strip() for seq in combined_seqs.split(',') if seq.strip()]
                            sequences_to_delete.update(sequences)
    
    return sequences_to_delete

def match_truncated_ids(original_fasta, truncated_ids):
    """
    Matches truncated sequence IDs to full headers in the original FASTA file.
    
    Args:
        original_fasta (str): Path to the original FASTA file.
        truncated_ids (set): Set of truncated sequence IDs to delete.
    
    Returns:
        set: Set of full headers corresponding to truncated IDs.
    """
    sequences_to_delete = set()
    headers_found = set()
    print("Matching truncated sequence IDs to full headers in the original FASTA...")
    
    try:
        for record in SeqIO.parse(original_fasta, "fasta"):
            header = record.description
            for truncated_id in truncated_ids:
                if truncated_id in header:
                    sequences_to_delete.add(header)
                    headers_found.add(truncated_id)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    
    # Identify truncated IDs not found
    missing_ids = truncated_ids - headers_found
    if missing_ids:
        print("\nWarning: The following sequence IDs were selected for deletion but were not found in the FASTA file:")
        for seq_id in missing_ids:
            print(seq_id)
    
    print(f"Total sequences to delete: {len(sequences_to_delete)}\n")
    return sequences_to_delete

def delete_sequences(original_fasta, sequences_to_delete, keep_fasta, deleted_fasta):
    """
    Deletes specified sequences from the original FASTA and writes two new FASTA files.
    
    Args:
        original_fasta (str): Path to the original FASTA file.
        sequences_to_delete (set): Set of full headers to delete.
        keep_fasta (str): Path to save the filtered FASTA.
        deleted_fasta (str): Path to save the deleted sequences FASTA.
    """
    to_keep = []
    to_delete = []
    total_sequences = 0

    print("Processing FASTA file for deletion...")
    try:
        for record in SeqIO.parse(original_fasta, "fasta"):
            total_sequences += 1
            if record.description in sequences_to_delete:
                to_delete.append(record)
            else:
                to_keep.append(record)
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        sys.exit(1)
    
    # Write the remaining sequences
    try:
        with open(keep_fasta, "w") as handle:
            SeqIO.write(to_keep, handle, "fasta")
        print(f"\nFiltered FASTA with remaining sequences saved to '{keep_fasta}'.")
    except Exception as e:
        print(f"Error writing kept FASTA file: {e}")
        sys.exit(1)
    
    # Write the deleted sequences
    try:
        with open(deleted_fasta, "w") as handle:
            SeqIO.write(to_delete, handle, "fasta")
        print(f"Deleted sequences FASTA saved to '{deleted_fasta}'.")
    except Exception as e:
        print(f"Error writing deleted FASTA file: {e}")
        sys.exit(1)
    
    print(f"\nTotal sequences processed: {total_sequences}")
    print(f"Sequences kept: {len(to_keep)}")
    print(f"Sequences deleted: {len(to_delete)}")

def main():
    args = parse_arguments()

    # Step 1: Extract gap-causing truncated sequence IDs from the report
    print("Extracting gap-causing sequences from the report...")
    gap_causers = extract_gap_causers(args.report)
    if not gap_causers:
        print("No gap-causing sequences found in the report. Exiting.")
        sys.exit(0)
    print(f"Total unique gap-causing truncated IDs found: {len(gap_causers)}\n")

    # Step 2: Ask user to select gap positions
    # First, gather all unique gap positions from the report
    try:
        excel_file = pd.ExcelFile(args.report)
    except Exception as e:
        print(f"Error reading Excel report: {e}")
        sys.exit(1)
    
    all_gap_positions = set()
    for sheet_name in excel_file.sheet_names:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
        all_gap_positions.update(df['Position'].dropna().astype(int).tolist())
    
    selected_positions = get_user_selected_positions(all_gap_positions)
    
    # Step 3: Collect sequences to delete based on selected positions
    sequences_to_delete_truncated = collect_sequences_to_delete(selected_positions, args.report)
    if not sequences_to_delete_truncated:
        print("No sequences to delete based on selected positions. Exiting.")
        sys.exit(0)
    print(f"Total unique truncated sequence IDs to delete: {len(sequences_to_delete_truncated)}\n")
    
    # Step 4: Match truncated IDs to full headers in the original FASTA
    sequences_to_delete_full = match_truncated_ids(args.fasta, sequences_to_delete_truncated)
    if not sequences_to_delete_full:
        print("No matching sequences found for deletion. Exiting.")
        sys.exit(0)
    
    # Step 5: Delete sequences and write output FASTA files
    delete_sequences(args.fasta, sequences_to_delete_full, args.keep_fasta, args.deleted_fasta)
    
    print("\nOperation completed successfully.")

if __name__ == "__main__":
    main()

