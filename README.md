Have you ever aligned a large set of sequences, only to find uninvited gaps in your multiple sequence alignment (MSA)?
Meet the MSA Gap Analyzer Scripts—a pair of Python tools designed to help you identify and remove problematic sequences causing gaps in your MSA.

gap_identifier.py
Takes an MSA file as input.
Pinpoints “gap causers” (the sequences responsible for introducing gaps).
Outputs a detailed report of gap positions and the specific sequences causing them.

delete_gap_causers.py
Lets you selectively remove identified gap-causing sequences from the original FASTA (used to generate the MSA).
Use these scripts to refine your alignment workflow and eliminate the hassle of unwanted gaps.
