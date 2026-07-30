from opengene_primers import (
    Primer,
    design_primers,
    find_primer_binding_sites,
    gc_content,
    melting_temperature,
    reverse_complement,
    simulate_pcr,
)


def test_sequence_metrics_and_primer_extension() -> None:
    primer = Primer("Tail", "ATGCGT", "GGCC")
    assert primer.full_sequence == "GGCCATGCGT"
    assert gc_content("GGAT") == 50.0
    assert melting_temperature("ATGC") == 12.0
    assert reverse_complement("ATGCRY") == "RYGCAT"


def test_binding_sites_include_both_orientations() -> None:
    sites = find_primer_binding_sites("AAATGCGCAT", "ATGC")
    assert [(site.start, site.end, site.strand) for site in sites] == [(2, 6, 1), (6, 10, -1)]


def test_pcr_simulation_incorporates_extensions() -> None:
    template = "AAATGCGTACGTTAGCGGTT"
    forward = Primer("F", "ATGCG", "GG")
    reverse = Primer("R", "AACCG", "TT")
    simulation = simulate_pcr(template, forward, reverse)
    assert len(simulation.products) == 1
    assert simulation.products[0].sequence == "GGATGCGTACGTTAGCGGTTAA"
    assert simulation.products[0].length == 22


def test_circular_pcr_crosses_origin() -> None:
    template = "TTAAACCCGGGATGCG"
    forward = Primer("F", "ATGCG")
    reverse = Primer("R", "TTAAA")
    reverse = Primer("R", "TTTAA")
    simulation = simulate_pcr(template, forward, reverse, circular=True)
    assert simulation.products
    assert simulation.products[0].sequence.startswith("ATGCG")


def test_simple_designer_returns_ranked_linear_pairs() -> None:
    template = "ATGCGTACGTTAGCGGTTACGATCGTACGATGCTAGCTAGCGTACGATCGTAGCTAGCATCGATCGTACGATCGTACG"
    pairs = design_primers(
        template,
        target_start=20,
        target_end=55,
        product_size_range=(50, 100),
        primer_length_range=(10, 14),
        tm_range=(24, 50),
        gc_range=(20, 80),
        search_window=20,
    )
    assert pairs
    assert all(50 <= pair.product_size <= 100 for pair in pairs)
    assert list(pairs) == sorted(pairs, key=lambda pair: (pair.score, pair.product_size, pair.forward.name, pair.reverse.name))
