#!/usr/bin/env python
"""
Stream IS.csv directly from GitHub (no local save) and build a synonym
lookup dictionary for IS elements.

Two dictionaries are produced:
  name_to_synonyms : canonical IS name -> list of synonym names
  synonym_to_name  : synonym name -> canonical IS name (for normalising
                      any IS name you found in your data back to its
                      canonical form)
"""

import csv
import io
import requests

RAW_URL = "https://raw.githubusercontent.com/thanhleviet/ISfinder-sequences/master/IS.csv"


def load_is_table(url: str = RAW_URL) -> list[list[str]]:
    """Fetch the CSV as text and parse it in memory. Nothing touches disk."""
    resp = requests.get(url)
    resp.raise_for_status()
    buffer = io.StringIO(resp.text)
    reader = csv.reader(buffer)
    return list(reader)


class UnionFind:
    def __init__(self):
        self.parent = {}

    def find(self, x):
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        rx, ry = self.find(x), self.find(y)
        if rx != ry:
            self.parent[rx] = ry


def build_synonym_dicts(rows: list[list[str]]):
    uf = UnionFind()

    for row in rows:
        if len(row) < 6:
            continue

        is_name = row[1].strip()
        synonym_field = row[4].strip()
        related_is = row[5].strip()

        if not is_name:
            continue

        uf.find(is_name)  # ensure it's registered even with no links

        if synonym_field and synonym_field.upper() != "NA":
            for syn in synonym_field.split(","):
                syn = syn.strip()
                if syn:
                    uf.union(is_name, syn)

        if related_is and related_is.upper() != "NA":
            uf.union(is_name, related_is)

    # group every name by its connected component
    groups = {}
    for name in uf.parent:
        root = uf.find(name)
        groups.setdefault(root, set()).add(name)

    name_to_synonyms = {}
    synonym_to_name = {}

    for group in groups.values():
        if len(group) < 2:
            continue
        # pick a stable canonical name (shortest, then alphabetical)
        canonical = sorted(group, key=lambda n: (len(n), n))[0]
        others = sorted(group - {canonical})
        name_to_synonyms[canonical] = others
        for member in group:
            if member != canonical:
                synonym_to_name[member] = canonical

    return name_to_synonyms, synonym_to_name

def get_synonyms(is_name: str, name_to_synonyms: dict, synonym_to_name: dict) -> list[str]:
    """
    Return all synonym names for a given IS name (excluding is_name itself).
    Works whether is_name is the canonical name or one of its synonyms.
    """
    canonical = synonym_to_name.get(is_name, is_name)
    syns = set(name_to_synonyms.get(canonical, []))
    if canonical != is_name:
        syns.add(canonical)
    syns.discard(is_name)
    return sorted(syns)


def create_dict(is_names: list[str], name_to_synonyms: dict, synonym_to_name: dict) -> dict[str, list[str]]:
    """
    Given a list of IS names, return a dict mapping each name to its list of
    synonyms. IS names with no real synonyms (nothing but themselves) are
    excluded from the result.
    """
    result = {}
    for name in is_names:
        syns = get_synonyms(name, name_to_synonyms, synonym_to_name)
        if syns:
            result[name] = syns
    return result


def final_dict(is_names: list[str], url: str = RAW_URL) -> dict[str, list[str]]:
    """
    All-in-one: fetch the ISfinder table, build the synonym lookup, and
    return a dict mapping each name in is_names to its list of synonyms.
    IS names with no synonyms are excluded from the result.

    This is the only function you need to call from your notebook/qmd:

        import IS_dict as isd
        is_dictionary = isd.get_canonical_dict(["IS26", "IS1358A", "IS150"])
    """
    rows = load_is_table(url)
    name_to_synonyms, synonym_to_name = build_synonym_dicts(rows)
    return create_dict(is_names, name_to_synonyms, synonym_to_name)


if __name__ == "__main__":
    # Simplest usage: just call final_dict with your list of IS names.
    # IS names with no real synonyms are excluded from the result.
    your_found_is_names = ["IS26", "IS15DI", "IS1358A", "ISEc53"]
    is_dictionary = final_dict(your_found_is_names)
    print(is_dictionary)