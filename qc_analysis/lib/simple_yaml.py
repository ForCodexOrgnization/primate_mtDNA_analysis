"""Restricted, standard-library-only YAML reader for repository QC configuration.

This is deliberately *not* a full YAML implementation.  It supports only the
mapping/scalar subset used by ``config/qc_preprocessing.yaml``: indentation
based mappings, empty mappings, simple scalars, quoted strings, and comments
outside quoted strings.  Use PyYAML for YAML documents needing sequences,
anchors, multiline values, or other YAML features.
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple


def parse_scalar(value: str) -> object:
    """Convert a scalar from the limited configuration subset."""
    value = value.strip()
    if value == "" or value.lower() in {"null", "none", "~"}:
        return None
    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def strip_yaml_comment(line: str) -> str:
    """Remove a comment marker unless it appears in a quoted string."""
    in_single = in_double = False
    for index, char in enumerate(line):
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single:
            in_double = not in_double
        elif char == "#" and not in_single and not in_double:
            return line[:index]
    return line


def read_simple_yaml(path: Path) -> Dict[str, object]:
    """Read the restricted mapping-only YAML subset used by QC configuration."""
    path = Path(path)
    root: Dict[str, object] = {}
    stack: List[Tuple[int, Dict[str, object]]] = [(-1, root)]
    implicit_mappings: set[int] = set()
    previous_indent = -1
    previous_was_mapping = True
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            if "\t" in raw_line[: len(raw_line) - len(raw_line.lstrip())]:
                raise ValueError(f"Unsupported tab indentation in {path}:{line_number}")
            line = strip_yaml_comment(raw_line.rstrip("\n"))
            if not line.strip():
                continue
            indent = len(line) - len(line.lstrip(" "))
            text = line.strip()
            if ":" not in text:
                raise ValueError(f"Unsupported YAML line in {path}:{line_number}: {raw_line.rstrip()}")
            if indent > previous_indent and not previous_was_mapping:
                raise ValueError(f"Invalid YAML indentation in {path}:{line_number}")
            while stack and indent <= stack[-1][0]:
                stack.pop()
            if not stack:
                raise ValueError(f"Invalid YAML indentation in {path}:{line_number}")
            key, value = text.split(":", 1)
            key = key.strip().strip('"\'')
            if not key:
                raise ValueError(f"Empty YAML key in {path}:{line_number}")
            value = value.strip()
            parent = stack[-1][1]
            if value == "":
                child: Dict[str, object] = {}
                parent[key] = child
                implicit_mappings.add(id(child))
                stack.append((indent, child))
                previous_was_mapping = True
            elif value == "{}":
                parent[key] = {}
                previous_was_mapping = False
            else:
                parent[key] = parse_scalar(value)
                previous_was_mapping = False
            previous_indent = indent

    def normalize_empty_implicit(mapping: Dict[str, object]) -> None:
        for key, value in mapping.items():
            if isinstance(value, dict):
                if id(value) in implicit_mappings and not value:
                    mapping[key] = None
                else:
                    normalize_empty_implicit(value)

    normalize_empty_implicit(root)
    return root
