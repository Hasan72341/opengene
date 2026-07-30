from opengene_restriction import (
    RestrictionEnzyme,
    Topology,
    digest,
    reverse_complement,
    scan_restriction_sites,
    unique_cutters,
)


def test_iupac_reverse_complement_and_site_scanning() -> None:
    enzyme = RestrictionEnzyme("AmbI", "GCNGC", 2, 3)
    assert reverse_complement("ARYN") == "NRYT"
    sites = scan_restriction_sites("TTGCA GCTT", [enzyme])
    assert len(sites) == 1
    assert (sites[0].start, sites[0].end, sites[0].strand) == (2, 7, 1)
    assert (sites[0].cut_forward, sites[0].cut_reverse) == (4, 5)


def test_non_palindromic_site_transforms_cuts_on_reverse_strand() -> None:
    enzyme = RestrictionEnzyme("TypeIIS", "GGTCTC", 7, 11)
    sites = scan_restriction_sites("GAGACCA", [enzyme])
    assert len(sites) == 1
    site = sites[0]
    assert site.strand == -1
    assert (site.cut_forward, site.cut_reverse) == (-5, -1)


def test_linear_digest_has_expected_fragments_and_unique_cutter() -> None:
    enzyme = RestrictionEnzyme("EcoRI", "GAATTC", 1, 5)
    result = digest("AAAAGAATTCCCC", [enzyme])
    assert result.cut_positions == (5,)
    assert [fragment.length for fragment in result.fragments] == [5, 8]
    assert [fragment.sequence for fragment in result.fragments] == ["AAAAG", "AATTCCCC"]
    assert unique_cutters("AAAAGAATTCCCC", [enzyme]) == (enzyme,)


def test_circular_digest_detects_origin_spanning_site_and_rotates_fragment() -> None:
    enzyme = RestrictionEnzyme("EcoRI", "GAATTC", 1, 5)
    sequence = "TTCGGGGAA"
    sites = scan_restriction_sites(sequence, [enzyme], topology=Topology.CIRCULAR)
    assert len(sites) == 1
    assert sites[0].start == 6
    assert sites[0].wraps_origin is True
    result = digest(sequence, [enzyme], topology="circular")
    assert result.cut_positions == (7,)
    assert len(result.fragments) == 1
    assert result.fragments[0].length == len(sequence)
    assert result.fragments[0].sequence == sequence[7:] + sequence[:7]


def test_linear_type_iis_outside_cut_is_reported_without_invalid_fragment() -> None:
    enzyme = RestrictionEnzyme("BsaI", "GGTCTC", 7, 11)
    result = digest("GGTCTC", [enzyme])
    assert len(result.sites) == 1
    assert result.cut_positions == ()
    assert [fragment.sequence for fragment in result.fragments] == ["GGTCTC"]
    assert result.warnings
