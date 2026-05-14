"""
msaclean — Clean and curate multiple sequence alignments for phylogenetics.

Quick start
-----------
>>> from msaclean import clean_alignment, clean_alignment_file, iterative_clean
>>> from Bio import AlignIO
>>>
>>> # One-call file cleaning
>>> report = clean_alignment_file(
...     "GH7_aligned.fasta", "GH7_clean.fasta",
...     whitelist={"BAB64565.3", "AAY83390.3"},
...     max_gap_fraction=0.95,
... )
>>> print(report)
>>>
>>> # Object-level cleaning
>>> aln = AlignIO.read("GH7_aligned.fasta", "fasta")
>>> clean, report = clean_alignment(aln, whitelist={"BAB64565.3"})
>>>
>>> # Iterative realignment + curation (requires MAFFT)
>>> result = iterative_clean(
...     "GH7_raw.fasta", output_fasta="GH7_curated.fasta",
...     whitelist={"BAB64565.3", "AAY83390.3"},
...     stop_at_n=20, remove_per_round=2,
... )
>>> print(result.log_table)
"""

from .core import (
    clean_alignment,
    clean_alignment_file,
    filter_empty_seqs,
    filter_gappy_seqs,
    protect_sequences,
    gappiest_sequences,
)
from .iterative import iterative_clean, IterativeCleanResult, IterationResult

__version__ = "0.1.0"
__all__ = [
    "clean_alignment",
    "clean_alignment_file",
    "filter_empty_seqs",
    "filter_gappy_seqs",
    "protect_sequences",
    "gappiest_sequences",
    "iterative_clean",
    "IterativeCleanResult",
    "IterationResult",
]
