#!/bin/python

import pytest
from eris.pipeline import Locus 

class FakeGenome:
    """Minimal indexable stand-in for GenomeAssembly."""
    def __init__(self, contigs: dict[str, str]):
        self._contigs = contigs

    def __getitem__(self, contig_id):
        return self._contigs[contig_id]


def test_single_contig_locus_slices_correctly():
    genome = FakeGenome({"ctg1": "AAAACCCCGGGGTTTT"})
    locus = Locus(
        id="l1", contig="ctg1", start=4, end=8,
        targets=[], passengers=[], upstream_flanks=[], downstream_flanks=[],
    )
    assert locus.extract_sequence(genome) == "CCCC"


def test_multi_contig_locus_should_not_return_whole_contigs():
    """
    Regression test for extract_sequence: a stitched locus must return only
    the region within [start, end] on each fragment, not the full contig.
    """
    genome = FakeGenome({
        "ctgA": "AAAAAAAAAA" + "CCCC",   # 10bp junk + 4bp real locus
        "ctgB": "GGGG" + "TTTTTTTTTT",   # 4bp real locus + 10bp junk
    })
    locus = Locus(
        id="l2", contig="ctgA|ctgB", start=10, end=4,
        targets=[], passengers=[], upstream_flanks=[], downstream_flanks=[],
    )
    seq = locus.extract_sequence(genome)
    assert len(seq) == 8, (
        f"got {len(seq)}bp back, expected 8bp — extract_sequence is "
        "returning full contig sequences instead of slicing each fragment"
    )
    assert seq == "CCCCGGGG"