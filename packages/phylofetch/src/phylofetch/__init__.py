"""
phylofetch — Fetch CDS sequences and taxonomy lineages from NCBI.

Quick start
-----------
>>> from phylofetch import fetch_cds, fetch_taxonomy, batch_fetch_cds
>>>
>>> # Single accession
>>> seq = fetch_cds("BAB64565.3")
>>> print(seq[:60])
>>>
>>> # Batch download to FASTA
>>> accs = ["BAB64565.3", "AAY83390.3", "XP_006965674.1"]
>>> batch_fetch_cds(accs, "my_sequences.fasta")
>>>
>>> # Taxonomy
>>> tax = fetch_taxonomy("BAB64565.3")
>>> print(tax["phylum"], tax["class"])
"""

from .cds import fetch_cds, batch_fetch_cds
from .taxonomy import fetch_taxonomy, batch_fetch_taxonomy
from .utils import set_email, RateLimiter

__version__ = "0.1.0"
__all__ = [
    "fetch_cds",
    "batch_fetch_cds",
    "fetch_taxonomy",
    "batch_fetch_taxonomy",
    "set_email",
    "RateLimiter",
]
