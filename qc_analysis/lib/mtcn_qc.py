"""Cohort QC mtCN: 2 * mean mtDNA depth / median nuclear depth.

The zero-nuclear-coverage behavior retains the historical HOM QC convention:
when nuclear median depth is exactly zero, use 1x as the denominator.
Missing or negative coverage is not rescued and yields an unavailable mtCN.
"""
from __future__ import annotations

import math


def finite_number(raw):
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def mtcn_from_coverage(mt_mean, nuclear_median):
    """Return (mtCN, adjusted_nuclear_depth, nuclear_zero_adjusted)."""
    mt = finite_number(mt_mean)
    nuclear = finite_number(nuclear_median)
    if mt is None or nuclear is None or mt < 0 or nuclear < 0:
        return None, None, False
    adjusted = 1.0 if nuclear == 0 else nuclear
    return 2.0 * mt / adjusted, adjusted, nuclear == 0


def mtcn_from_summary(row):
    """Read standardized collection-summary fields, not historical mtcn_median."""
    return mtcn_from_coverage(
        row.get('mt_mean_coverage'),
        row.get('nuclear_median_coverage'),
    )
