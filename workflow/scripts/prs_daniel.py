#!/usr/bin/env python3
"""PRS calculation -- Python 3.9 port of "PRS313_calculation/prs.Daniel.R".

Original R script: Daniel Eriksson 2022.

This script calculates a PRS for an individual. In short, it takes a VCF file
and a score file, multiplies the genotypes with the corresponding weights
to sum a total score.

This script assumes one subject per gVCF file, and uses only a single score
file. This allows for missing calls to be identified, which in turn enables
imputation of missing scores and individual sensitivity analysis.

The PRS calculation has a few obstacles that makes the calculation (and the
script), a bit more complicated. The transition from genotyping to sequencing-
based PRS has necessitated the introduction of proxies and multi-allelic sites.
This, in turn, requires processing of alternative ways to encode the same
genetic alteration. To handle synonymous terminology, an individual's two
alleles are handled separately and the score file offers multiple ways that
alleles can yield scores. Before calculation of the grand total, any form of
double counting is actively avoided. Furthermore, the gVCF format contains
blocks of reference calls that overlap score file positions and must not be
interpreted as missing values. gVCF blocks are searched for overlaps with the
score file.


Input files, each named outright on the command line:
1. G-zipped gVCF file with extracted SNPs/variants (see NOTE on gzip below)
2. Polygenic risk score info file = score file

The marker range file (snps_ranges.tdf.<build>) is NOT an input here. It belongs
to the extraction step in vcf_pipeline.sh, which passes it to "bcftools view
-R"; by the time this script runs, the VCF has already been subset to those
positions. Note that the range file and the score file must be of the same
genome build, and that this script cannot check that for you -- a mismatch
simply looks like a sample with almost everything missing.

Output files, each written to the exact path the caller names:
1. --detailed: genotypes and scores for each position
2. --total: total PRS and number of missing calls:
   SumPRS        - An individuals total PRS, defined as the sum of all
                   genotype * weight products
   CountNA       - Number of positions/markers missing data
   MaxNeg        - The lowest possible score the individual could get if all
                   missing values were protective alleles.
   MaxPos        - The largest possible score the individual could get if all
                   missing values were risk alleles.
   MaxRange      - The difference between MaxPos and MaxNeg, corresponds to the
                   total weight of scores not assigned due to missing data.
   ExpectedScore - The PRS the individual would get if all missing scores are
                   replaced by the expected score, given the marker weights and
                   the allele frequency in the general population.


Usage:

    python3 prs_daniel.py -s <score_file> -v <vcf_file> \
                          -d <detailed_csv> -t <total_csv>

      -s, --score-file FILE   the score file
      -v, --vcf FILE          the VCF to score
      -d, --detailed FILE     path of the per-marker output CSV
      -t, --total FILE        path of the summary output CSV

All four are required and may be given in any order. Every path names one
specific file: no directory is searched for an input, and no name is derived
for an output. Missing parent directories of an output are created; an existing
file at that path is overwritten. Nothing about the run -- who started it,
when, which sample it is -- is known to this script: the calling pipeline owns
all of that, including what the two output files are called.

Requires Python 3.9+ and the standard library only (no vcfR equivalent needed;
the VCF parsing this script needs is limited to the fixed columns, the INFO END
key and the GT field of the first sample).


NOTES ON DELIBERATE DIFFERENCES FROM prs.Daniel.R
-------------------------------------------------
This is otherwise a faithful port, including the quirks documented inline below.
The changes are these -- 1 and 2 because the R script's behaviour cannot work
or cannot be reproduced; 3 and 4 because the R behaviour is wrong. Only 4 can
change a score, and --r-compat restores the R behaviour of 3 and 4:

0. The run identity -- username, date and project name -- is gone, along with
   the "<name>results_<user>_<date>" folder and the "<name>_<kind>_scores_<user>
   _<date>.csv" file names the R script built from it. This script builds no
   names at all: both output files are given as full paths, because the
   pipeline that calls it tracks the run. --r-compat does NOT bring the old
   naming back.
1. Both input files are named outright. The R script instead rebuilds the VCF
   name from the run identity and globs "score_file.*csv" out of a directory --
   so with one score file per genome build in it, which one gets used depends on
   alphabetical order. Here, picking the file is the caller's job, and a
   directory where a file is expected is an error rather than a search.
2. gzip is detected from the file's magic bytes, not its name. Depending on the
   bcftools version, "bcftools annotate -O v --output ...recode.vcf.gz" writes
   either an uncompressed VCF under a .gz name (the -O v wins) or a real bgzip
   file (the extension wins). Both are read.
3. The symbolic allele <NON_REF> is reported in full in the VCF_Alt column of
   detailed_scores.csv. prs.Daniel.R packs the allele and the reference into one
   "<allele>_<ref>" string and splits it again on "_", which truncates
   <NON_REF> to <NON. Scores are identical either way, because no effect allele
   contains an underscore. Pass --r-compat to reproduce the truncation and get
   byte-identical output.
4. A genotype counts as missing if ANY of its alleles is a no-call, whatever the
   spelling. prs.Daniel.R matches a fixed list of five strings ("./.", "./1",
   "1/.", "0/.", "0|.") before it normalises "|" to "/", so every phased
   spelling and every allele index above 1 escapes it -- and an escaped
   genotype is not just kept, it is miscounted: R's sort() drops the missing
   allele index, so ".|1" is read as the homozygote "1/1" and scores TWO effect
   alleles at a marker where one allele is unknown. Here such a marker is
   reported as missing instead, so it is counted in CountNA and carried through
   the sensitivity analysis. This is the one difference that moves scores, and
   it moves them only for VCFs that contain partially called genotypes in a
   spelling the R list misses. A position whose only record is dropped this way
   is also protected from the gVCF block expansion, so that an overlapping
   block cannot turn a no-call back into a confident reference call.

Output formatting follows R's write.csv(row.names = FALSE) -- quoted header and
character fields, bare numbers, bare NA, 15 significant digits -- so results can
be diffed directly against the R script's output. With --r-compat the two are
byte-identical; see tests/compare_r_python.sh.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import sys
from typing import Dict, List, Optional, Sequence, Tuple

# Settings -----------------------------------------------------------

# Default is True. Set to False only for detailed studies of proxy concordance.
AVOID_PROXY_DOUBLE_COUNTING = True

# The literal list from prs.Daniel.R. It is applied there BEFORE "|" is
# normalised to "/", so it catches only unphased spellings and only allele
# index 1: ".|1", "1|.", ".|.", "./2", "2/." and a bare "." all slip through.
# Worse, a genotype that slips through is not merely kept -- R's sort() drops
# the missing allele index, so ".|1" is read as the homozygote "1/1" and scores
# two effect alleles where the data says one allele is unknown.
# Used only under --r-compat.
R_INCOMPLETE_GENOTYPES = ("./.", "./1", "1/.", "0/.", "0|.")

# Separators between the alleles of a genotype: "/" unphased, "|" phased.
_GT_SEPARATOR = re.compile(r"[/|]")


def is_incomplete_genotype(genotype: Optional[str], r_compat: bool = False) -> bool:
    """True when a genotype is absent, or any one of its alleles is a no-call.

    A call is only usable if every allele of it is known. Anything else is
    missing data for that marker, which is what lets the sensitivity analysis
    account for it. This covers every spelling: unphased and phased, any
    ploidy, and any allele index -- "./.", ".|.", "./1", ".|1", "1|.", "0/.",
    "./0", "2|.", a bare "." and an empty field.
    """
    if genotype is None:
        return True
    genotype = genotype.strip()
    if genotype in ("", "."):
        return True
    if r_compat:
        return genotype in R_INCOMPLETE_GENOTYPES
    return any(allele.strip() in ("", ".") for allele in _GT_SEPARATOR.split(genotype))


# Column layout of the score file. The names in the file's header row are
# ignored; columns are taken by position, exactly like read.csv(col.names = ...).
SCORE_COLUMNS: Sequence[Tuple[str, str]] = (
    ("Chr", "int"),
    ("Pos", "int"),
    ("Ref", "str"),
    ("Alt", "str"),
    ("RSID", "str"),
    ("Effect_Allele", "str"),
    ("Weight", "float"),
    ("Effect_allele_frequency", "float"),
    ("MarkerStatus", "str"),
)

DETAILED_COLUMNS = [name for name, _ in SCORE_COLUMNS] + [
    "ID",
    "ChrPos",
    "VCF_Ref",
    "VCF_Alt",
    "Number_of_effect_alleles",
    "Product",
]

SUMMARY_COLUMNS = [
    "SumPRS",
    "CountNA",
    "MaxNeg",
    "MaxPos",
    "MaxRange",
    "ExpectedScore",
    "DiffSumVsExpected",
]

Row = Dict[str, object]


def log(message: str) -> None:
    print(message, flush=True)


# Reading helpers -----------------------------------------------------------


def _na(value: str) -> Optional[str]:
    """Mimic read.csv's default na.strings = "NA" on an unquoted field."""
    value = value.strip()
    if value == "" or value == "NA":
        return None
    return value


def _as_int(value: str) -> Optional[int]:
    """as.numeric() on a field that holds a whole number; NA on failure."""
    text = _na(value)
    if text is None:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return int(number)


def _as_float(value: str) -> Optional[float]:
    text = _na(value)
    if text is None:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _split_csv_line(line: str) -> List[str]:
    """Split one CSV line, honouring double quotes the way read.csv does."""
    fields: List[str] = []
    field: List[str] = []
    in_quotes = False
    index = 0
    while index < len(line):
        char = line[index]
        if in_quotes:
            if char == '"':
                if index + 1 < len(line) and line[index + 1] == '"':
                    field.append('"')
                    index += 1
                else:
                    in_quotes = False
            else:
                field.append(char)
        elif char == '"':
            in_quotes = True
        elif char == ",":
            fields.append("".join(field))
            field = []
        else:
            field.append(char)
        index += 1
    fields.append("".join(field))
    return fields


# Read and prepare Score file -----------------------------------------------


def prepare_output(path: str, what: str) -> str:
    """Make an output file writable: reject a directory, create the parent."""
    if os.path.isdir(path):
        raise SystemExit("%s must be a file, not a directory: %s" % (what, path))
    parent = os.path.dirname(os.path.abspath(path))
    try:
        os.makedirs(parent, exist_ok=True)
    except OSError as error:
        raise SystemExit("Cannot create the directory for %s: %s" % (what, error))
    return path


def require_file(path: str, what: str) -> str:
    """Every input is a named file. Say so clearly when one is not."""
    if os.path.isdir(path):
        raise SystemExit("%s must be a file, not a directory: %s" % (what, path))
    if not os.path.isfile(path):
        raise SystemExit("%s not found: %s" % (what, path))
    return path


def read_score_file(score_file: str) -> List[Row]:
    log("Reading in score info file.")
    require_file(score_file, "Score file")

    rows: List[Row] = []
    with open(score_file, "r", encoding="utf-8-sig", newline="") as handle:
        header = handle.readline()
        if not header:
            raise SystemExit("Score file %s is empty." % score_file)
        for line_number, line in enumerate(handle, start=2):
            line = line.rstrip("\r\n")
            if not line:
                continue
            fields = _split_csv_line(line)
            if len(fields) != len(SCORE_COLUMNS):
                raise SystemExit(
                    "Score file %s line %d: expected %d columns, found %d. "
                    "Expected layout: %s"
                    % (
                        score_file,
                        line_number,
                        len(SCORE_COLUMNS),
                        len(fields),
                        ", ".join(name for name, _ in SCORE_COLUMNS),
                    )
                )
            row: Row = {}
            for (name, kind), raw in zip(SCORE_COLUMNS, fields):
                if kind == "int":
                    row[name] = _as_int(raw)
                elif kind == "float":
                    row[name] = _as_float(raw)
                else:
                    row[name] = _na(raw)
            # Give IDs to the alleles and the positions
            row["ID"] = "%s:%s:%s:%s" % (
                _r_chr(row["Chr"]),
                _r_chr(row["Pos"]),
                _r_chr(row["Ref"]),
                _r_chr(row["Alt"]),
            )
            row["ChrPos"] = "%s:%s" % (_r_chr(row["Chr"]), _r_chr(row["Pos"]))
            rows.append(row)

    if not rows:
        raise SystemExit("Score file %s holds no markers." % score_file)
    log("Score file: %s (%d markers)" % (score_file, len(rows)))
    return rows


def _r_chr(value: object) -> str:
    """as.character() as used by paste(): NA becomes the string "NA"."""
    if value is None:
        return "NA"
    if isinstance(value, float):
        return format_number(value)
    return str(value)


# Read and prepare VCF file -------------------------------------------------


def _open_vcf(path: str):
    """Open a VCF that may or may not actually be gzipped, whatever it is named."""
    with open(path, "rb") as probe:
        magic = probe.read(2)
    if magic == b"\x1f\x8b":
        import gzip

        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


_END_INFO = re.compile(r"(?:^|;)END=([0-9]+)")


def read_vcf(path: str, r_compat: bool = False) -> Tuple[List[Row], set]:
    log("Reading in VCF file.")

    rows: List[Row] = []
    with _open_vcf(path) as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            line = line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                continue

            chrom = _as_int(re.sub(r"^chr", "", fields[0]))
            pos = _as_int(fields[1])
            info = fields[7]

            # Add info about END.
            # gVCFs contain blocks spanning many positions; record the end of
            # each block so that the block can be matched against the score file.
            # (The R script pulled END out with a chain of greedy gsub()s; this
            # reads the INFO key itself, which is the same thing for a gVCF but
            # is not fooled by other keys whose name ends in "END".)
            end_match = _END_INFO.search(info)
            end = int(end_match.group(1)) if end_match else pos

            # Add genotypes: GT of the first sample, as extract.gt(element='GT').
            genotype: Optional[str] = None
            if len(fields) > 9 and fields[8]:
                keys = fields[8].split(":")
                if "GT" in keys:
                    values = fields[9].split(":")
                    index = keys.index("GT")
                    if index < len(values):
                        genotype = _na(values[index])

            rows.append(
                {
                    "CHROM": chrom,
                    "POS": pos,
                    "ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "END": end if end is not None else pos,
                    "Genotype": genotype,
                }
            )

    log("Read %d variant records." % len(rows))
    if not rows:
        raise SystemExit(
            "No variant records in %s -- the extraction step probably failed. "
            "Check the size of that file; vcf_pipeline.sh does not check "
            "whether bcftools succeeded." % path
        )

    # Remove variant calls with missing or incomplete genotypes
    def chrpos(row: Row) -> str:
        return "%s:%s" % (_r_chr(row["CHROM"]), _r_chr(row["POS"]))

    before = len(rows)
    no_call = {chrpos(row) for row in rows
               if is_incomplete_genotype(row["Genotype"], r_compat)}
    rows = [row for row in rows if not is_incomplete_genotype(row["Genotype"], r_compat)]
    if before != len(rows):
        log("Dropped %d records with a missing or incomplete genotype."
            % (before - len(rows)))
    for row in rows:
        row["Genotype"] = str(row["Genotype"]).replace("|", "/", 1)
        row["BlockLength"] = (
            row["END"] - row["POS"]
            if row["END"] is not None and row["POS"] is not None
            else None
        )

    # Remove duplicates, leaving one
    rows = _unique_by(rows, lambda row: (row["CHROM"], row["POS"]))
    log("Kept %d variant records with usable genotypes." % len(rows))
    # Positions left with no call at all. They must stay missing: an overlapping
    # gVCF block would otherwise turn a no-call into a confident reference call,
    # which is the opposite of what the data says. Positions that also carry a
    # usable record are not protected -- that record is the call.
    return rows, no_call - {chrpos(row) for row in rows}


def _unique_by(rows: List[Row], key) -> List[Row]:
    """duplicated() semantics: keep the first row of every key."""
    seen = set()
    kept = []
    for row in rows:
        marker = key(row)
        if marker in seen:
            continue
        seen.add(marker)
        kept.append(row)
    return kept


# Convert gVCF blocks into localized reference calls -------------------------


def expand_gvcf_blocks(
    vcf_rows: List[Row], score_rows: List[Row], no_call: Optional[set] = None
) -> List[Row]:
    """Turn every score position covered by a gVCF block into a reference call.

    For every position in the score file that is covered by a gVCF block, add an
    additional row to the VCF describing this particular position as a reference
    call. Without this, confidently-called reference positions inside a block
    would be counted as missing data.
    """
    log("Expanding gVCF blocks.")
    no_call = no_call or set()

    by_chrom: Dict[Optional[int], List[Row]] = {}
    for row in vcf_rows:
        by_chrom.setdefault(row["CHROM"], []).append(row)

    def overlapping_row(chrom: Optional[int], pos: Optional[int]) -> Optional[Row]:
        # Does the score file position overlap a block in the VCF? If so, return
        # that block's row; if not, return None.
        if chrom is None or pos is None:
            return None
        for row in by_chrom.get(chrom, ()):
            if row["END"] is not None and row["POS"] is not None:
                if row["END"] >= pos and row["POS"] <= pos:
                    return row
        return None

    expanded = list(vcf_rows)
    added = 0
    protected = 0
    for score_row in score_rows:
        chrom = score_row["Chr"]
        pos = score_row["Pos"]
        if score_row["ChrPos"] in no_call:
            protected += 1
            continue
        block = overlapping_row(chrom, pos)
        if block is None:
            continue

        new_row = dict(block)
        # Adjust the position to match the overlapping position in the score file
        new_row["POS"] = pos
        # R for reference, but not specified which base because not every ref
        # base is given in gVCF blocks, only the ref of the first base.
        new_row["REF"] = "R"
        new_row["ALT"] = "<NON_REF>"
        new_row["ID"] = "%s : %s :R:NA" % (_r_chr(new_row["CHROM"]), _r_chr(pos))
        expanded.append(new_row)
        by_chrom.setdefault(new_row["CHROM"], []).append(new_row)
        added += 1

    for row in expanded:
        row["ChrPos"] = "%s:%s" % (_r_chr(row["CHROM"]), _r_chr(row["POS"]))

    # Real variant rows come first, so they win over a synthesised reference call
    # at the same position.
    expanded = _unique_by(expanded, lambda row: row["ChrPos"])
    log("Added %d block-derived reference calls." % added)
    if protected:
        log("Left %d no-call position(s) missing rather than reading an "
            "overlapping block as a reference call." % protected)
    return expanded


# Convert number genotypes into allele name genotypes ------------------------


def add_letter_genotypes(vcf_rows: List[Row]) -> None:
    """For example: 1/3 --> AA/AATAT, given REF A and ALT AA,AAT,AATAT."""
    for row in vcf_rows:
        # Temporarily fuse all alleles including the reference:
        # "REF,ALT1,ALT2,...,ALTn"; REF is item number 1.
        alleles = ("%s,%s" % (row["REF"], row["ALT"])).split(",")

        indices: List[int] = []
        for part in str(row["Genotype"]).split("/"):
            try:
                indices.append(int(float(part)))
            except ValueError:
                # as.numeric(".") is NA, and R's sort() drops NAs -- so a
                # half-called genotype collapses to the called allele alone.
                continue
        indices.sort()

        letters = []
        for index in indices:
            # +1 in R; here the reference is simply element 0.
            letters.append(alleles[index] if 0 <= index < len(alleles) else "NA")
        letter_genotype = "/".join(letters)

        # Replace the marker of a spanning deletion
        row["LetterGenotype"] = letter_genotype.replace("*", "", 1)


# Split genotypes to enable simplification of variant descriptions ------------


def simplify_incl_ref(allele: str, ref: str, r_compat: bool = False) -> Tuple[str, str]:
    """Simplify a variant description, e.g. Ref CTT Alt CT -> Ref CT Alt C."""
    if r_compat:
        # prs.Daniel.R packs the allele and the reference into a single
        # "<allele>_<ref>" string and splits it again on "_". An allele that
        # contains an underscore itself -- in practice only the symbolic
        # <NON_REF> -- is therefore truncated to "<NON". This is score
        # neutral, because no effect allele contains an underscore, but it is
        # visible in the VCF_Alt column of detailed_scores.csv. Enable it only
        # to reproduce the R output byte for byte.
        packed = "%s_%s" % (allele, ref)
        allele = packed.split("_", 1)[0]
        ref = packed.rsplit("_", 1)[-1]
    while len(allele) > 1 and len(ref) > 1 and allele[-1] == ref[-1]:
        # As long as the last base is the same, remove it from both
        allele = allele[:-1]
        ref = ref[:-1]
    return allele, ref


def build_sides(vcf_rows: List[Row], r_compat: bool = False) -> Tuple[List[Row], List[Row]]:
    """Split the VCF into two "monoploid" tables, one per allele of the subject.

    Correct scoring of the variants relies on them being described in the most
    simple way, in order for the variants to be comparable to the score file.
    Both alleles are simplified independently with regards to the reference,
    which means that from this point a single position (Chr:Pos) may have two
    different VCF-Ref alleles in the left and the right table.
    """
    log("Simplify variant descriptions")

    # Sort the table
    ordered = sorted(
        vcf_rows,
        key=lambda row: (
            row["CHROM"] is None,
            row["CHROM"] if row["CHROM"] is not None else 0,
            row["POS"] is None,
            row["POS"] if row["POS"] is not None else 0,
        ),
    )

    left: List[Row] = []
    right: List[Row] = []
    for marker_id, row in enumerate(ordered, start=1):
        letter_genotype = str(row["LetterGenotype"])
        # Split every genotype A/B in two, one table counting the risk allele(s)
        # in the left part and one in the right part of the genotype.
        allele_a = letter_genotype.split("/", 1)[0]
        allele_b = letter_genotype.rsplit("/", 1)[-1]

        for table, allele in ((left, allele_a), (right, allele_b)):
            simplified_allele, simplified_ref = simplify_incl_ref(
                allele, str(row["REF"]), r_compat
            )
            table.append(
                {
                    "MarkerID": marker_id,
                    "ChrPos": row["ChrPos"],
                    "REF": simplified_ref,
                    "Allele_X": simplified_allele,
                }
            )
    return left, right


# Fuse VCF file with score file to join genotypes and weights ----------------


def build_subject_scores(
    score_rows: List[Row], left: List[Row], right: List[Row]
) -> List[Row]:
    left_by_pos = {row["ChrPos"]: row for row in left}
    right_by_pos = {row["ChrPos"]: row for row in right}

    def side_table(vcf_by_pos: Dict[str, Row]) -> List[Row]:
        table: List[Row] = []
        for score_row in score_rows:
            row = dict(score_row)
            # Insert the REF and ALT from the VCF into the new table
            vcf_row = vcf_by_pos.get(str(score_row["ChrPos"]))
            if vcf_row is None:
                row["VCF_Ref"] = None
                row["VCF_Alt"] = None
                # No call at this position: the marker is missing data (NA),
                # which is different from a call that simply does not match.
                row["t1"] = None
            else:
                row["VCF_Ref"] = vcf_row["REF"]
                row["VCF_Alt"] = vcf_row["Allele_X"]
                # Is the VCF variant (ALT) describing an effect variant listed
                # in the score file?
                row["t1"] = int(
                    row["VCF_Alt"] == score_row["Effect_Allele"]
                    and row["VCF_Ref"] == score_row["Ref"]
                )
            table.append(row)

        # Sort by number of risk alleles (NA last), so that the effect allele is
        # the row that survives deduplication of a position.
        table.sort(key=lambda row: (row["t1"] is None, -(row["t1"] or 0)))
        table = _unique_by(table, lambda row: row["ChrPos"])
        # Sort the table by ChrPos again
        table.sort(
            key=lambda row: (
                row["Chr"] is None,
                row["Chr"] if row["Chr"] is not None else 0,
                row["Pos"] if row["Pos"] is not None else 0,
            )
        )
        return table

    left_scores = side_table(left_by_pos)
    right_scores = side_table(right_by_pos)

    # Checkpoint: order and content in both parallel data sets are the same
    assert [row["ChrPos"] for row in left_scores] == [
        row["ChrPos"] for row in right_scores
    ], "The two monoploid score tables went out of sync."

    # Fuse score sheets into one
    subject_scores: List[Row] = []
    for left_row, right_row in zip(left_scores, right_scores):
        row = dict(left_row)
        del row["t1"]
        if left_row["t1"] is None or right_row["t1"] is None:
            row["Number_of_effect_alleles"] = None
        else:
            row["Number_of_effect_alleles"] = left_row["t1"] + right_row["t1"]
        subject_scores.append(row)
    return subject_scores


# Avoid double counting of multi-allelic sites and proxies -------------------


def _feed_back(rows: List[Row], totals: Dict[str, Optional[int]]) -> None:
    """Feed the number of effect alleles per site back to the main score sheet."""
    for row in rows:
        status = row["MarkerStatus"]
        if status in totals:
            row["Number_of_effect_alleles"] = totals[status]


def collapse_multiallelic(rows: List[Row]) -> List[Row]:
    """Sum the effect alleles of each multi-allelic site, then keep one row.

    A site is missing as a whole if any of its rows is missing, and a site can
    never contribute more than the two alleles the subject has.
    """
    groups: Dict[str, List[Optional[int]]] = {}
    for row in rows:
        status = row["MarkerStatus"]
        if isinstance(status, str) and "MULTIALLELIC" in status:
            groups.setdefault(status, []).append(row["Number_of_effect_alleles"])

    totals: Dict[str, Optional[int]] = {}
    for status, counts in groups.items():
        if any(count is None for count in counts):
            totals[status] = None
        else:
            totals[status] = min(2, sum(count for count in counts if count is not None))
    _feed_back(rows, totals)

    # Remove duplicates, leaving one, for each multi-allelic site, to avoid
    # double counting the weight
    seen = set()
    kept: List[Row] = []
    for row in rows:
        status = row["MarkerStatus"]
        multiallelic = isinstance(status, str) and "MULTIALLELIC" in status
        if multiallelic and status in seen:
            continue
        seen.add(status)
        kept.append(row)
    return kept


def collapse_proxies(rows: List[Row]) -> List[Row]:
    """Pair each replaced site with its proxy and keep one row of the pair."""
    groups: Dict[str, List[Optional[int]]] = {}
    for row in rows:
        status = row["MarkerStatus"]
        if isinstance(status, str) and "PROXY" in status:
            # Give proxies and the sites they replace the same ID, forming pairs
            key = status.replace("REPLACED_BY_PROXY", "PROXY")
            groups.setdefault(key, []).append(row["Number_of_effect_alleles"])

    totals: Dict[str, Optional[int]] = {}
    for key, counts in groups.items():
        called = [count for count in counts if count is not None]
        # If both the replaced position and the proxy are missing, the pair is
        # missing; otherwise take the call that is available.
        totals[key] = min(2, max(called)) if called else None
    _feed_back(rows, totals)

    # Remove the replaced sites, leaving the proxies, to avoid double counting
    return [
        row
        for row in rows
        if not (
            isinstance(row["MarkerStatus"], str)
            and "REPLACED_BY_PROXY" in row["MarkerStatus"]
        )
    ]


# Calculate PRS and run the individual sensitivity analysis ------------------


def calculate(rows: List[Row]) -> Row:
    log("Calculating summary statistics.")

    for row in rows:
        count = row["Number_of_effect_alleles"]
        weight = row["Weight"]
        row["Product"] = None if count is None or weight is None else count * weight

    # The score disregarding any NA's
    score = math.fsum(
        row["Product"] for row in rows if isinstance(row["Product"], float)
    )
    missing = [row for row in rows if row["Number_of_effect_alleles"] is None]

    summary: Row = {"SumPRS": score, "CountNA": len(missing)}

    log("Performing sensitivity analysis.")

    # Estimate what the score would be if an individual's missing calls were all
    # 1) risk alleles, or 2) protective alleles.
    if missing:
        weights = [row["Weight"] for row in missing if row["Weight"] is not None]
        # The sum of all positive scores is the worst case scenario
        max_pos = math.fsum(w for w in weights if w > 0) * 2 + score
        max_neg = math.fsum(w for w in weights if w < 0) * 2 + score
    else:
        # Kept from prs.Daniel.R: with nothing missing these are reported as 0,
        # not as the score itself.
        max_pos = 0.0
        max_neg = 0.0

    summary["MaxNeg"] = max_neg
    summary["MaxPos"] = max_pos
    summary["MaxRange"] = max_pos - max_neg

    # Estimate what the score would be if an individual's missing calls were all
    # the expected value, given the effect allele frequency in the population.
    if missing:
        if any(
            row["Effect_allele_frequency"] is None or row["Weight"] is None
            for row in missing
        ):
            summary["ExpectedScore"] = None
        else:
            # 2*p*q + 2*p^2 is the expected number of effect alleles (= 2p).
            # Summed in one go, as R's sum() does, so that the result agrees
            # with the R script to the last digit.
            products = []
            for row in missing:
                p = row["Effect_allele_frequency"]
                q = 1 - p
                products.append((2 * p * q + 2 * p ** 2) * row["Weight"])
            summary["ExpectedScore"] = math.fsum(products) + score
    else:
        summary["ExpectedScore"] = score

    summary["DiffSumVsExpected"] = (
        None
        if summary["ExpectedScore"] is None
        else summary["SumPRS"] - summary["ExpectedScore"]
    )
    return summary


# Write result files ---------------------------------------------------------


def format_number(value: float) -> str:
    """as.character() on a double: 15 significant digits, no trailing ".0"."""
    if math.isnan(value):
        return "NaN"
    if math.isinf(value):
        return "Inf" if value > 0 else "-Inf"
    if value == 0:
        # R prints a bare 0; Python would render a negative zero as "-0".
        return "0"
    return "%.15g" % value


def _format_field(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, str):
        return '"%s"' % value.replace('"', '""')
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return format_number(value)
    return str(value)


def write_r_csv(path: str, columns: Sequence[str], rows: Sequence[Row]) -> None:
    """write.csv(..., row.names = FALSE): quoted header and strings, bare NA."""
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write(",".join('"%s"' % column for column in columns) + "\n")
        for row in rows:
            handle.write(
                ",".join(_format_field(row.get(column)) for column in columns) + "\n"
            )


# Entry point ----------------------------------------------------------------


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate a polygenic risk score for one subject from a gVCF and a "
            "score file (Python port of prs.Daniel.R). Writes the two result "
            "CSVs to the paths given by --detailed and --total."
        )
    )
    parser.add_argument(
        "-d",
        "--detailed",
        required=True,
        metavar="FILE",
        help="path of the per-marker output CSV to write",
    )
    parser.add_argument(
        "-t",
        "--total",
        required=True,
        metavar="FILE",
        help="path of the summary output CSV to write",
    )
    parser.add_argument(
        "-s",
        "--score-file",
        required=True,
        metavar="FILE",
        help="the score file: a 9-column CSV, one row per marker",
    )
    parser.add_argument(
        "-v",
        "--vcf",
        dest="vcf_file",
        required=True,
        metavar="FILE",
        help="the extracted VCF for one subject (.vcf or .vcf.gz)",
    )
    parser.add_argument(
        "--no-proxy-dedup",
        action="store_true",
        help=(
            "keep both the replaced site and its proxy; for detailed studies of "
            "proxy concordance only, as it double counts their weights"
        ),
    )
    parser.add_argument(
        "--r-compat",
        action="store_true",
        help=(
            "reproduce prs.Daniel.R byte for byte: its truncation of the "
            "symbolic allele <NON_REF> to <NON in the VCF_Alt column, and its "
            "fixed five-string list of missing genotypes, which misreads a "
            "partially called phased genotype such as '.|1' as a homozygote"
        ),
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    log("Start of Python script %s" % os.path.basename(__file__))

    if os.path.abspath(args.detailed) == os.path.abspath(args.total):
        raise SystemExit("--detailed and --total must name different files.")

    score_rows = read_score_file(args.score_file)

    vcf_path = require_file(args.vcf_file, "VCF file")
    log("VCF file: %s" % vcf_path)
    vcf_rows, no_call = read_vcf(vcf_path, args.r_compat)
    if args.r_compat:
        # prs.Daniel.R lets a block override a dropped call; keep that here.
        no_call = set()
    if not vcf_rows:
        # Legitimate, if extreme: every marker is simply missing data.
        log("WARNING: no usable genotypes -- every marker will be reported as "
            "missing.")

    vcf_rows = expand_gvcf_blocks(vcf_rows, score_rows, no_call)
    add_letter_genotypes(vcf_rows)
    left, right = build_sides(vcf_rows, args.r_compat)

    subject_scores = build_subject_scores(score_rows, left, right)
    subject_scores = collapse_multiallelic(subject_scores)
    if AVOID_PROXY_DOUBLE_COUNTING and not args.no_proxy_dedup:
        subject_scores = collapse_proxies(subject_scores)

    summary = calculate(subject_scores)

    log("Writing results files.")
    detailed = prepare_output(args.detailed, "--detailed")
    totals = prepare_output(args.total, "--total")
    write_r_csv(detailed, DETAILED_COLUMNS, subject_scores)
    write_r_csv(totals, SUMMARY_COLUMNS, [summary])

    log("SumPRS = %s (CountNA = %d)" % (format_number(summary["SumPRS"]), summary["CountNA"]))
    log("Wrote %s" % detailed)
    log("Wrote %s" % totals)
    return 0


if __name__ == "__main__":
    sys.exit(main())
