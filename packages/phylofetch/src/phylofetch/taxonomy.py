"""Fetch full taxonomy lineages from NCBI for protein/nucleotide accessions."""

import logging
from typing import Optional, Sequence, Union
from pathlib import Path

import pandas as pd
from Bio import Entrez

from .utils import _default_limiter, _ensure_email, with_retry

logger = logging.getLogger(__name__)

# Ranks we extract into the returned dict
_TARGET_RANKS = ("superkingdom", "kingdom", "phylum", "class",
                 "order", "family", "genus", "species")


def fetch_taxonomy(accession: str,
                   extra_ranks: Sequence[str] = (),
                   retries: int = 3) -> dict:
    """Fetch the full taxonomy lineage for a protein or nucleotide accession.

    Uses NCBI Entrez to retrieve the GenBank record, extract the TaxID,
    then resolves the lineage to common taxonomic ranks.

    Parameters
    ----------
    accession : str
        NCBI accession (protein or nucleotide).
    extra_ranks : list of str, optional
        Additional NCBI rank names to include in the output beyond the
        default set (superkingdom → species).
    retries : int
        Retry attempts on network failure. Default 3.

    Returns
    -------
    dict
        Keys: ``"accession"``, ``"taxid"``, ``"lineage_string"``, plus one key
        per rank in `_TARGET_RANKS` (``"superkingdom"``, ``"phylum"``, …).
        Values are ``None`` when a rank is absent for that organism.

    Examples
    --------
    >>> from phylofetch import set_email, fetch_taxonomy
    >>> set_email("you@example.com")
    >>> tax = fetch_taxonomy("BAB64565.3")
    >>> print(tax["phylum"], "/", tax["class"])
    'Metamonada / Parabasalia'
    """
    _ensure_email()
    _default_limiter.wait()

    all_ranks = list(_TARGET_RANKS) + list(extra_ranks)
    result = {"accession": accession, "taxid": None, "lineage_string": None}
    result.update({r: None for r in all_ranks})

    @with_retry
    def _fetch(acc: str) -> dict:
        # Get GenBank record to find the TaxID
        handle = Entrez.efetch(db="protein", id=acc,
                               rettype="gb", retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        taxid = None
        lin_str = None
        for rec in records:
            if "GBSeq_taxonomy" in rec:
                lin_str = rec["GBSeq_taxonomy"]
            # TaxID is in the features
            for feat in rec.get("GBSeq_feature-table", []):
                for qual in feat.get("GBFeature_quals", []):
                    if qual.get("GBQualifier_name") == "db_xref":
                        val = qual.get("GBQualifier_value", "")
                        if val.startswith("taxon:"):
                            taxid = val.split(":")[1]
                            break
                if taxid:
                    break

        if not taxid:
            logger.warning("No TaxID found for %s", acc)
            return result

        result["taxid"] = taxid
        result["lineage_string"] = lin_str

        # Resolve lineage via Taxonomy DB
        _default_limiter.wait()
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        tax_records = Entrez.read(handle)
        handle.close()

        if not tax_records:
            return result

        lineage = tax_records[0].get("LineageEx", [])
        for node in lineage:
            rank = node.get("Rank", "").lower()
            name = node.get("ScientificName", "")
            if rank in all_ranks:
                result[rank] = name

        # species from the record itself
        result["species"] = tax_records[0].get("ScientificName")

        return result

    return _fetch(accession)


def batch_fetch_taxonomy(accessions: Sequence[str],
                         output_excel: Optional[Union[str, Path]] = None,
                         output_tsv: Optional[Union[str, Path]] = None,
                         extra_ranks: Sequence[str] = (),
                         retries: int = 3,
                         verbose: bool = True) -> pd.DataFrame:
    """Fetch taxonomy for a list of accessions and return a tidy DataFrame.

    Parameters
    ----------
    accessions : list of str
        NCBI accession numbers.
    output_excel : str or Path, optional
        Save results to an Excel (.xlsx) file.
    output_tsv : str or Path, optional
        Save results to a tab-separated file.
    extra_ranks : list of str, optional
        Additional ranks to include.
    retries : int
        Per-accession retry attempts. Default 3.
    verbose : bool
        Print progress. Default ``True``.

    Returns
    -------
    pandas.DataFrame
        One row per accession with columns: accession, taxid,
        lineage_string, superkingdom, phylum, class, order, family,
        genus, species (and any extras).

    Examples
    --------
    >>> from phylofetch import set_email, batch_fetch_taxonomy
    >>> set_email("you@example.com")
    >>> accs = ["BAB64565.3", "XP_006965674.1", "WP_165634028.1"]
    >>> df = batch_fetch_taxonomy(accs, output_excel="taxonomy.xlsx")
    >>> df[["accession", "phylum", "class"]]
    """
    _ensure_email()
    rows = []

    for i, acc in enumerate(accessions, 1):
        if verbose:
            print(f"[{i}/{len(accessions)}] {acc} ...", end=" ", flush=True)
        try:
            row = fetch_taxonomy(acc, extra_ranks=extra_ranks, retries=retries)
            rows.append(row)
            if verbose:
                print(row.get("phylum") or "?")
        except Exception as exc:
            logger.error("Failed for %s: %s", acc, exc)
            rows.append({"accession": acc})
            if verbose:
                print(f"ERROR: {exc}")

    df = pd.DataFrame(rows)

    if output_excel:
        df.to_excel(output_excel, index=False)
        if verbose:
            print(f"\nSaved → {output_excel}")

    if output_tsv:
        df.to_csv(output_tsv, sep="\t", index=False)
        if verbose:
            print(f"Saved → {output_tsv}")

    return df
