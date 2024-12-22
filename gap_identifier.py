import argparse
import sys
from Bio import AlignIO
import pandas as pd
import math

def parse_arguments():
    parser = argparse.ArgumentParser(description="Identify gap causers (sequences without gaps) in an MSA and save the report as an Excel file with separate handling for large numbers of gap-causers.")
    parser.add_argument("-i", "--input", required=True, help="Input MSA file (e.g., FASTA, Clustal, Stockholm).")
    parser.add_argument("-f", "--format", required=False, default="fasta",
                        help="Format of the MSA file (default: fasta). Options include clustal, stockholm, etc.")
    parser.add_argument("-o", "--output", required=False, default="gap_causers_report.xlsx",
                        help="Output Excel file to save the gap causers report (default: gap_causers_report.xlsx).")
    parser.add_argument("-c", "--cutoff", type=int, required=False, default=1,
                        help="Minimum number of sequences without a gap to report per position (default: 1).")
    parser.add_argument("--small_chunk_size", type=int, required=False, default=1000,
                        help="Number of positions per sheet for small gap causers (default: 1000).")
    parser.add_argument("--large_chunk_size", type=int, required=False, default=1000,
                        help="Number of positions per sheet for large gap causers (default: 1000).")
    return parser.parse_args()

def identify_gap_causers(alignment, cutoff=1, threshold=2000):
    small_gap_causers = []
    large_gap_causers = []
    alignment_length = alignment.get_alignment_length()

    # Define gap characters
    gap_characters = ['-', '.', '?']  # Adjust if your MSA uses different gap symbols

    # Iterate through each position
    for pos in range(alignment_length):
        causer_sequences = [record.id for record in alignment if record.seq[pos] not in gap_characters]
        num_causers = len(causer_sequences)
        if num_causers >= cutoff:
            if num_causers <= threshold:
                small_gap_causers.append({
                    'Position': pos + 1,  # 1-based index
                    'Number_of_Gap_Causers': num_causers,
                    'Gap_Causers_Sequences': causer_sequences  # Keep as list
                })
            else:
                # To prevent excessively long strings, you might consider limiting the number of sequences concatenated
                # or handle them in a different manner. Here, we'll proceed as per your request.
                large_gap_causers.append({
                    'Position': pos + 1,  # 1-based index
                    'Number_of_Gap_Causers': num_causers,
                    'Gap_Causers_Sequences_Combined': ", ".join(causer_sequences)
                })
    return small_gap_causers, large_gap_causers

def expand_small_gap_causers(df_small, max_sequences=2000):
    if df_small.empty:
        return df_small  # Return empty DataFrame if no data

    # Determine the maximum number of gap causers in small positions
    max_causers = df_small['Gap_Causers_Sequences'].apply(len).max()

    # Ensure max_causers does not exceed the threshold
    max_causers = min(max_causers, max_sequences)

    # Create new columns for each gap causer
    for i in range(max_causers):
        column_name = f'Gap_Causers_Sequence_{i+1}'
        df_small[column_name] = df_small['Gap_Causers_Sequences'].apply(
            lambda x: x[i] if i < len(x) else ""
        )

    # Drop the original list column
    df_small = df_small.drop(columns=['Gap_Causers_Sequences'])
    return df_small

def write_excel_report(small_gap_causers, large_gap_causers, output_file, small_chunk_size=1000, large_chunk_size=1000, max_sequences=2000):
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        # Handle Small Gap Causers
        if small_gap_causers:
            df_small = pd.DataFrame(small_gap_causers)
            df_small = expand_small_gap_causers(df_small, max_sequences=max_sequences)

            # Determine the number of chunks/sheets needed
            total_small_positions = len(df_small)
            num_small_chunks = math.ceil(total_small_positions / small_chunk_size)

            for i in range(num_small_chunks):
                start = i * small_chunk_size
                end = start + small_chunk_size
                chunk_df = df_small.iloc[start:end].reset_index(drop=True)
                sheet_name = f'Small_Gap_Causers_{i+1}'

                # Check Excel's sheet name length limit (31 characters)
                if len(sheet_name) > 31:
                    sheet_name = f'Small_Gap_{i+1}'

                chunk_df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"Written sheet: {sheet_name} (Positions {start +1} to {min(end, total_small_positions)})")
        else:
            print("No small gap causers found meeting the specified criteria.")

        # Handle Large Gap Causers
        if large_gap_causers:
            df_large = pd.DataFrame(large_gap_causers)

            # Calculate the number of chunks/sheets needed
            total_large_positions = len(df_large)
            num_large_chunks = math.ceil(total_large_positions / large_chunk_size)

            for i in range(num_large_chunks):
                start = i * large_chunk_size
                end = start + large_chunk_size
                chunk_df = df_large.iloc[start:end].reset_index(drop=True)
                sheet_name = f'Large_Gap_Causers_{i+1}'

                # Check Excel's sheet name length limit (31 characters)
                if len(sheet_name) > 31:
                    sheet_name = f'Large_Gap_{i+1}'

                chunk_df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"Written sheet: {sheet_name} (Positions {start +1} to {min(end, total_large_positions)})")
        else:
            print("No large gap causers found meeting the specified criteria.")

    print(f"Gap causers report successfully saved to {output_file}")

def main():
    args = parse_arguments()

    try:
        alignment = AlignIO.read(args.input, args.format)
    except Exception as e:
        print(f"Error reading alignment file: {e}")
        sys.exit(1)

    print("Identifying gap causers...")
    small_gap_causers, large_gap_causers = identify_gap_causers(alignment, cutoff=args.cutoff, threshold=2000)

    print("Writing to Excel...")
    write_excel_report(
        small_gap_causers,
        large_gap_causers,
        args.output,
        small_chunk_size=args.small_chunk_size,
        large_chunk_size=args.large_chunk_size,
        max_sequences=2000
    )

if __name__ == "__main__":
    main()

