from pathlib import Path

import pytest

from kffractbias.analysis import AnalysisConfig, calculate_fractionation_bias


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_plot_outputs_are_created_and_stale_plots_are_removed(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    target = write(tmp_path / "target.bed", "chr2\t0\t10\tt1\nchr10\t0\t10\tt2\n")
    query = write(tmp_path / "query.bed", "chrA\t0\t10\tq1\n")
    anchors = write(tmp_path / "pairs.anchors", "t1\tq1\t10\nt2\tq1\t10\n")
    output_dir = tmp_path / "output"
    plotted = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=output_dir,
            window_size=1,
        )
    )

    assert plotted.pdf_path is not None and plotted.pdf_path.stat().st_size > 0
    assert plotted.png_path is not None and plotted.png_path.stat().st_size > 0

    unplotted = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=output_dir,
            window_size=1,
            make_plot=False,
        )
    )
    assert unplotted.pdf_path is None
    assert unplotted.png_path is None
    assert not (output_dir / "kffractbias.plot.pdf").exists()
    assert not (output_dir / "kffractbias.plot.png").exists()


def test_plot_failure_preserves_previous_outputs(tmp_path: Path, monkeypatch) -> None:
    target = write(tmp_path / "target.bed", "chr1\t0\t10\tt1\n")
    query = write(tmp_path / "query.bed", "chrA\t0\t10\tq1\n")
    anchors = write(tmp_path / "pairs.anchors", "t1\tq1\t10\n")
    output_dir = tmp_path / "output"
    original = calculate_fractionation_bias(
        AnalysisConfig(
            synteny_path=anchors,
            synteny_format="jcvi",
            target_bed=target,
            query_bed=query,
            output_dir=output_dir,
            window_size=1,
            make_plot=False,
        )
    )
    previous_genes = original.genes_path.read_bytes()
    previous_summary = original.summary_path.read_bytes()

    def fail_plot(*_args, **_kwargs):
        raise RuntimeError("synthetic plot failure")

    monkeypatch.setattr("kffractbias.plotting.plot_windows", fail_plot)
    with pytest.raises(RuntimeError, match="synthetic plot failure"):
        calculate_fractionation_bias(
            AnalysisConfig(
                synteny_path=anchors,
                synteny_format="jcvi",
                target_bed=target,
                query_bed=query,
                output_dir=output_dir,
                window_size=1,
                make_plot=True,
            )
        )

    assert original.genes_path.read_bytes() == previous_genes
    assert original.summary_path.read_bytes() == previous_summary


def test_many_sequences_are_paginated_without_losing_panels(tmp_path):
    pytest.importorskip("matplotlib")
    import re

    import matplotlib.pyplot as plt

    from kffractbias.io import Gene
    from kffractbias.plotting import plot_windows
    from kffractbias.profiles import RetentionProfile

    profile = RetentionProfile(
        {f"t{i}": [Gene(f"t{i}", 0, 3, f"gene{i}")] for i in range(8)},
        [f"q{i}" for i in range(13)],
        {},
        1,
        1,
    )
    pdf, png = tmp_path / "plot.pdf", tmp_path / "plot.png"
    metadata = plot_windows(
        profile.window_rows(), pdf, png, target_name="target", query_name="query", window_size=1
    )
    assert metadata["panels"] == 16
    assert metadata["pdf_pages"] == 3
    assert len(re.findall(rb"/Type /Page\b", pdf.read_bytes())) == 3
    assert png.stat().st_size > 0
    assert not plt.get_fignums()
