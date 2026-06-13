#!/usr/bin/env python3
"""
Compile a CSV of ciliary proteins from the Human Protein Atlas (HPA).

Background
----------
Hansen et al., Cell 2025 ("Intrinsic heterogeneity of primary cilia revealed
through spatial proteomics", S0092-8674(25)01029-3) expanded the HPA with
antibody-based subciliary annotations. They annotated four primary-cilium
regions and identified 715 cilium-localizing genes:
    "Primary cilium"                  ~378 genes
    "Primary cilium tip"              ~117 genes
    "Primary cilium transition zone"  ~100 genes
    "Basal body"                      ~438 genes
    (union ~715)

Those four labels are HPA subcellular-location categories, so we can pull the
current ciliary proteome straight out of the HPA query/download API and join on
all the other per-gene annotations HPA holds (coordinates, tissue expression,
disease, etc.).

What this script does
---------------------
1. Queries HPA's search_download API once per cilium compartment, filtering on
   subcellular location.
2. Unions the four result sets, de-duplicating by Ensembl gene ID and recording
   which compartment(s) each gene was found in.
3. Writes one tidy CSV.

It uses only the Python standard library (urllib, json, csv) so it runs anywhere
proteinatlas.org is reachable (laptop, GenomeDK login node, etc.) with no pip
install.

API reference: https://www.proteinatlas.org/about/help/dataaccess
Data version at time of writing: HPA v25.1 (Ensembl v109 / GRCh38).
"""

import csv
import json
import sys
import time
import urllib.parse
import urllib.request

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

API = "https://www.proteinatlas.org/api/search_download.php"

# The four HPA subcellular-location categories that define a "ciliary" protein.
# Edit this list to broaden/narrow (e.g. drop "Basal body", or add
# "Centriolar satellite" / "Centrosome" if you want the wider MTOC context).
CILIUM_LOCATIONS = [
    "Primary cilium",
    "Primary cilium tip",
    "Primary cilium transition zone",
    "Basal body",
    "Centrosome",            # wider MTOC context
    "Centriolar satellite",  # wider MTOC context
]

# Expected per-compartment counts from the paper, for a quick sanity check.
EXPECTED = {
    "Primary cilium": 378,
    "Primary cilium tip": 117,
    "Primary cilium transition zone": 100,
    "Basal body": 438,
}

# Columns to request, as (api_code, output_header).
# Full list of codes: https://www.proteinatlas.org/about/help/dataaccess
BASE_COLUMNS = [
    ("g",      "Gene"),
    ("gs",     "Gene synonym"),
    ("eg",     "Ensembl"),
    ("gd",     "Gene description"),
    ("up",     "Uniprot"),
    ("chr",    "Chromosome"),
    ("chrp",   "Position"),                 # chromosomal location (GRCh38)
    ("pc",     "Protein class"),
    ("upbp",   "Biological process"),
    ("up_mf",  "Molecular function"),
    ("di",     "Disease involvement"),
    ("pe",     "Evidence"),
    ("scl",    "Subcellular location"),     # all annotated locations
    ("scml",   "Subcellular main location"),
    ("scal",   "Subcellular additional location"),
    ("relce",  "Reliability (IF)"),         # immunofluorescence reliability
    ("rnats",  "RNA tissue specificity"),
    ("rnatd",  "RNA tissue distribution"),
    ("rnatss", "RNA tissue specificity score"),
    ("rnatsm", "RNA tissue specific nTPM"),  # names the top-expressed tissues
    ("rnascs", "RNA single cell type specificity"),
    ("rnascsm","RNA single cell type specific nCPM"),
    ("prcts",  "Protein cell type specificity"),  # DVP, incl. ciliated cells
]

# Set True to also append the full consensus per-tissue RNA matrix (51 tissues,
# nTPM). Gives the quantitative "expressed in which tissues" answer column by
# column. Set False for a leaner file (the rnatsm summary above still names the
# enriched tissues).
INCLUDE_ALL_TISSUE_NTPM = False

CONSENSUS_TISSUES = [
    "adipose_tissue", "adrenal_gland", "amygdala", "appendix", "basal_ganglia",
    "blood_vessel", "bone_marrow", "breast", "cerebellum", "cerebral_cortex",
    "cervix", "choroid_plexus", "colon", "duodenum", "endometrium_1",
    "epididymis", "esophagus", "fallopian_tube", "gallbladder", "heart_muscle",
    "hippocampal_formation", "hypothalamus", "kidney", "liver", "lung",
    "lymph_node", "midbrain", "ovary", "pancreas", "parathyroid_gland",
    "pituitary_gland", "placenta", "prostate", "rectum", "retina",
    "salivary_gland", "seminal_vesicle", "skeletal_muscle", "skin_1",
    "small_intestine", "smooth_muscle", "spinal_cord", "spleen", "stomach_1",
    "testis", "thymus", "thyroid_gland", "tongue", "tonsil", "urinary_bladder",
    "vagina",
]

OUTPUT_CSV = "ciliary_proteins_hpa.csv"

# --------------------------------------------------------------------------- #
# Build the full column list
# --------------------------------------------------------------------------- #

COLUMNS = list(BASE_COLUMNS)
if INCLUDE_ALL_TISSUE_NTPM:
    for t in CONSENSUS_TISSUES:
        COLUMNS.append((f"t_RNA_{t}", f"Tissue RNA {t} [nTPM]"))

COLUMN_CODES = ",".join(code for code, _ in COLUMNS)


# --------------------------------------------------------------------------- #
# Fetch
# --------------------------------------------------------------------------- #

def fetch_location(location, retries=3, pause=3.0):
    """Return a list of record dicts for one subcellular-location query."""
    params = {
        "search": f"subcell_location:{location}",
        "format": "json",
        "columns": COLUMN_CODES,
        "compress": "no",
    }
    url = f"{API}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(
        url, headers={"User-Agent": "ciliary-proteome-compiler/1.0"}
    )
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                raw = resp.read().decode("utf-8")
            return json.loads(raw)
        except Exception as e:  # noqa: BLE001
            last_err = e
            sys.stderr.write(
                f"  attempt {attempt}/{retries} for '{location}' failed: {e}\n"
            )
            time.sleep(pause)
    raise RuntimeError(f"Could not fetch '{location}': {last_err}")


def normalize(value):
    """Flatten HPA's list-valued JSON fields into a single string cell."""
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(str(v) for v in value)
    return str(value)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def main():
    headers = [hdr for _, hdr in COLUMNS]
    # The JSON keys returned by HPA are the human-readable headers above.
    merged = {}  # ensembl id -> row dict

    for loc in CILIUM_LOCATIONS:
        sys.stderr.write(f"Fetching '{loc}' ...\n")
        records = fetch_location(loc)
        n = len(records)
        exp = EXPECTED.get(loc)
        note = f" (paper reported ~{exp})" if exp else ""
        sys.stderr.write(f"  -> {n} genes{note}\n")
        if n == 0:
            sys.stderr.write(
                "  WARNING: 0 results. The 'subcell_location' field token or the\n"
                "  category spelling may have changed in a newer HPA release.\n"
            )
        for rec in records:
            ens = rec.get("Ensembl") or rec.get("Gene") or json.dumps(rec)
            row = merged.setdefault(ens, {})
            for _, hdr in COLUMNS:
                if hdr not in row or not row[hdr]:
                    row[hdr] = normalize(rec.get(hdr))
            comps = row.setdefault("_compartments", set())
            comps.add(loc)

    if not merged:
        sys.stderr.write("No data retrieved. Aborting.\n")
        sys.exit(1)

    out_headers = ["Cilium compartment(s)"] + headers
    rows = []
    for row in merged.values():
        comps = sorted(row.pop("_compartments", set()))
        row["Cilium compartment(s)"] = "; ".join(comps)
        rows.append(row)

    # Stable sort: by chromosome then gene name.
    def sort_key(r):
        return (str(r.get("Chromosome", "")), str(r.get("Gene", "")))

    rows.sort(key=sort_key)

    with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_headers, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    sys.stderr.write(
        f"\nWrote {len(rows)} unique ciliary genes to {OUTPUT_CSV}\n"
        f"(paper's union was ~715; small differences reflect HPA version drift).\n"
    )


if __name__ == "__main__":
    main()
