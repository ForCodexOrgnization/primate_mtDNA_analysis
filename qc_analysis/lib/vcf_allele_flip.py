"""Helpers for safe REF/ALT orientation changes in lifted VCF records.

These rules mirror the REF_ALT_FLIP logic used by the NUMT consensus workflow;
keep status names and conservative field handling in sync when either workflow is
updated.
"""

from __future__ import annotations

import re
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

BASES = {"A", "C", "G", "T"}
REF_MATCH = "REF_MATCH"
ALT_REF_FLIP = "ALT_REF_FLIP"
FLIP_UNSUPPORTED_VARIANT_TYPE = "FLIP_UNSUPPORTED_VARIANT_TYPE"
TARGET_REF_NOT_SOURCE_REF_OR_ALT = "TARGET_REF_NOT_SOURCE_REF_OR_ALT"
UNMAPPED = "UNMAPPED"


def is_simple_biallelic_snv(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and ref.upper() in BASES and alt.upper() in BASES and "," not in alt


def classify_ref_alt_relationship(source_ref: str, source_alt: str, target_ref: str, enable_ref_alt_flip: bool = True) -> str:
    if target_ref.upper() == source_ref.upper():
        return REF_MATCH
    if target_ref.upper() == source_alt.upper():
        if enable_ref_alt_flip and is_simple_biallelic_snv(source_ref, source_alt):
            return ALT_REF_FLIP
        return FLIP_UNSUPPORTED_VARIANT_TYPE
    return TARGET_REF_NOT_SOURCE_REF_OR_ALT


def split_meta(line: str) -> Dict[str, str]:
    m = re.match(r"##(?:INFO|FORMAT)=<(.+)>", line.strip())
    if not m:
        return {}
    text = m.group(1)
    out: Dict[str, str] = {}
    key = ""
    val = ""
    in_key = True
    in_quote = False
    pairs: List[Tuple[str, str]] = []
    for ch in text + ",":
        if in_key:
            if ch == "=":
                in_key = False
            else:
                key += ch
        else:
            if ch == '"':
                in_quote = not in_quote
                val += ch
            elif ch == "," and not in_quote:
                pairs.append((key.strip(), val.strip().strip('"')))
                key, val, in_key = "", "", True
            else:
                val += ch
    for k, v in pairs:
        out[k] = v
    return out


def swap_ref_alt_array(value: str) -> Tuple[str, bool]:
    vals = value.split(",")
    if len(vals) == 2:
        return f"{vals[1]},{vals[0]}", True
    return ".", False


def invert_af(value: str) -> Tuple[str, bool]:
    vals = value.split(",")
    if len(vals) != 1:
        return ".", False
    try:
        af = float(vals[0])
    except ValueError:
        return ".", False
    if not 0.0 <= af <= 1.0:
        return ".", False
    new = 1.0 - af
    return ("{:.6g}".format(new), True)


def af_from_ad(old_ad: Optional[str]) -> Optional[str]:
    if not old_ad:
        return None
    vals = old_ad.split(",")
    if len(vals) != 2:
        return None
    try:
        ref_depth = float(vals[0])
        alt_depth = float(vals[1])
    except ValueError:
        return None
    total = ref_depth + alt_depth
    if total <= 0:
        return None
    return "{:.6g}".format(ref_depth / total)


def flip_gt(gt: str) -> str:
    if gt in {"", "."}:
        return gt
    parts = re.split(r"([/|])", gt)
    flipped = []
    for token in parts:
        if token == "0":
            flipped.append("1")
        elif token == "1":
            flipped.append("0")
        else:
            flipped.append(token)
    return "".join(flipped)


def _flip_number_g(value: str, gt: Optional[str]) -> Tuple[str, bool]:
    vals = value.split(",")
    alleles = [] if not gt else [a for a in re.split(r"[/|]", gt) if a != "."]
    ploidy = len(alleles)
    if len(vals) == 2 and ploidy == 1:  # haploid: 0,1 -> 1,0
        return f"{vals[1]},{vals[0]}", True
    if len(vals) == 3 and ploidy in {0, 2}:  # diploid VCF order 00,01,11 -> 11,01,00
        return f"{vals[2]},{vals[1]},{vals[0]}", True
    return ".", False


def flip_sample_fields(format_keys: Sequence[str], sample_value: str, format_numbers: Mapping[str, str]) -> Tuple[str, List[str]]:
    vals = sample_value.split(":")
    while len(vals) < len(format_keys):
        vals.append("")
    original = dict(zip(format_keys, vals))
    cleared: List[str] = []
    old_ad = original.get("AD")
    for i, key in enumerate(format_keys):
        val = vals[i]
        if val in {"", "."}:
            continue
        number = format_numbers.get(key)
        if key == "GT":
            vals[i] = flip_gt(val)
        elif key in {"AD", "FAD", "F1R2", "F2R1"} or number == "R":
            vals[i], ok = swap_ref_alt_array(val)
            if not ok:
                cleared.append(key)
        elif key == "AF":
            derived = af_from_ad(old_ad)
            if derived is not None:
                vals[i] = derived
            else:
                vals[i], ok = invert_af(val)
                if not ok:
                    cleared.append(key)
        elif key == "SB":
            sb = val.split(",")
            if len(sb) == 4:
                vals[i] = ",".join(sb[2:4] + sb[0:2])
            else:
                vals[i] = "."
                cleared.append(key)
        elif number == "G" or key in {"PL", "GL", "GP"}:
            vals[i], ok = _flip_number_g(val, original.get("GT"))
            if not ok:
                cleared.append(key)
    return ":".join(vals[: len(format_keys)]), cleared


def parse_info(info: str) -> MutableMapping[str, Optional[str]]:
    d: MutableMapping[str, Optional[str]] = {}
    if not info or info == ".":
        return d
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
        else:
            d[item] = None
    return d


def format_info(info: Mapping[str, Optional[str]]) -> str:
    if not info:
        return "."
    parts = []
    for k, v in info.items():
        parts.append(k if v is None else f"{k}={v}")
    return ";".join(parts)


def transform_info_fields(info_str: str, info_numbers: Mapping[str, str]) -> Tuple[str, List[str]]:
    info = parse_info(info_str)
    dropped: List[str] = []
    for key in list(info.keys()):
        val = info[key]
        number = info_numbers.get(key)
        if val is None:
            continue
        if number == "R":
            info[key], ok = swap_ref_alt_array(val)
            if not ok:
                del info[key]
                dropped.append(key)
        elif number in {"A", "G"}:
            del info[key]
            dropped.append(key)
    if dropped:
        info["LIFTOVER_DROPPED_INFO_FIELDS"] = ",".join(dropped)
    return format_info(info), dropped
