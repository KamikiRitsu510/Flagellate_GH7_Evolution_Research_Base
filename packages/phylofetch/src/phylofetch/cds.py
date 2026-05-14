"""Fetch CDS (coding DNA sequences) from NCBI for protein accessions.

NCBI stores CDS in different places depending on the record type:
- Nucleotide records: sequence is the CDS itself
- Protein records: CDS coordinates are stored in a ``coded_by`` qualifier
  pointing to a nucleotide record

This module handles both cases transparently.
"""

import re
import logging
from pathlib import Path
from typing import Optional, Sequence, Union

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .utils import _default_limiter, _ensure_email, with_retry

logger = logging.getLogger(__name__)


# ── Core fetch ────────────────────────────────────────────────────────────────

def fetch_cds(accession: str,
              retries: int = 3) -> Optional[str]:
    """Fetch the CDS nucleotide sequence for a protein or nucleotide accession.

    The function first searches the nucleotide database. If nothing is found
    it falls back to the protein database and extracts the CDS coordinates
    from the ``coded_by`` qualifier.

    Parameters
    ----------
    accession : str
        An NCBI accession number (protein or nucleotide), e.g. ``"BAB64565.3"``.
    retries : int
        Number of retry attempts on network errors. Default 3.

    Returns
    -------
    str or None
        The CDS sequence as a plain string (ACGT), or ``None`` if not found.

    Examples
    --------
    >>> from phylofetch import set_email, fetch_cds
    >>> set_email("you@example.com")
    >>> seq = fetch_cds("BAB64565.3")
    >>> len(seq)
    1035
    """
    _ensure_email()
    _default_limiter.wait()

    @with_retry
    def _fetch(acc: str) -> Optional[str]:
        # ── Try nucleotide DB first ──
        handle = Entrez.esearch(db="nucleotide", term=acc)
        nt_result = Entrez.read(handle)
        handle.close()

        if nt_result["Count"] != "0":
            nuc_id = nt_result["IdList"][0]
            h = Entrez.efetch(db="nucleotide", id=nuc_id,
                              rettype="fasta", retmode="text")
            rec = SeqIO.read(h, "fasta")
            h.close()
            return str(rec.seq)

        # ── Fall back to protein DB, extract coded_by coordinates ──
        handle = Entrez.esearch(db="protein", term=acc)
        prot_result = Entrez.read(handle)
        handle.close()

        if prot_result["Count"] == "0":
            logger.warning("Accession not found in NCBI: %s", acc)
            return None

        prot_id = prot_result["IdList"][0]
        h = Entrez.efetch(db="protein", id=prot_id,
                           rettype="gb", retmode="text")
        prot_rec = SeqIO.read(h, "genbank")
        h.close()

        for feature in prot_rec.features:
            if feature.type == "CDS" and "coded_by" in feature.qualifiers:
                coded_by = feature.qualifiers["coded_by"][0]
                # Format: "NM_012345.1:100..500" or "complement(NM_012345.1:100..500)"
                match = re.search(r"([A-Z_0-9.]+):(\d+)\.\.(\d+)", coded_by)
                if match:
                    nuc_acc   = match.group(1)
                    seq_start = int(match.group(2))
                    seq_stop  = int(match.group(3))
                    _default_limiter.wait()
                    h2 = Entrez.efetch(db="nucleotide", id=nuc_acc,
                                       rettype="fasta", retmode="text",
                                       seq_start=seq_start, seq_stop=seq_stop)
                    nuc_rec = SeqIO.read(h2, "fasta")
                    h2.close()
                    return str(nuc_rec.seq)

        logger.warning("CDS not found in protein record for: %s", acc)
        return None

    return _fetch(accession)


# ── Batch fetch ───────────────────────────────────────────────────────────────

def batch_fetch_cds(accessions: Sequence[str],
                    output_fasta: Optional[Union[str, Path]] = None,
                    retries: int = 3,
                    verbose: bool = True) -> dict:
    """Fetch CDS sequences for a list of accessions.

    Parameters
    ----------
    accessions : list of str
        NCBI accession numbers to fetch.
    output_fasta : str or Path, optional
        If provided, write results to this FASTA file. Missing sequences are
        silently skipped.
    retries : int
        Per-accession retry attempts. Default 3.
    verbose : bool
        Print progress to stdout. Default ``True``.

    Returns
    -------
    dict
        ``{accession: sequence_string}`` for successful fetches.
        Accessions for which no sequence was found are absent from the dict.

    Examples
    --------
    >>> from phylofetch import set_email, batch_fetch_cds
    >>> set_email("you@example.com")
    >>> accs = ["BAB64565.3", "AAY83390.3", "BAB64553.1"]
    >>> seqs = batch_fetch_cds(accs, "flagellate_cds.fasta")
    >>> print(f"Downloaded {len(seqs)}/{len(accs)} sequences")
    """
    _ensure_email()
    results: dict = {}
    records = []

    for i, acc in enumerate(accessions, 1):
        if verbose:
            print(f"[{i}/{len(accessions)}] Fetching {acc} ...", end=" ", flush=True)
        seq = fetch_cds(acc, retries=retries)
        if seq:
            results[acc] = seq
            records.append(SeqRecord(Seq(seq), id=acc, description=""))
            if verbose:
                print(f"OK ({len(seq)} bp)")
        else:
            if verbose:
                print("NOT FOUND")

    if output_fasta and records:
        output_fasta = Path(output_fasta)
        with open(output_fasta, "w") as fh:
            SeqIO.write(records, fh, "fasta")
        if verbose:
            print(f"\nSaved {len(records)} sequences → {output_fasta}")

    return results
