from pathlib import Path

import pytest

from kffractbias.jcvi import run_pairwise_synteny, run_self_synteny


@pytest.fixture(autouse=True)
def no_external_tools(monkeypatch):
    monkeypatch.setattr("kffractbias.jcvi.preflight_tools", lambda _aligner: None)


def write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_run_pairwise_synteny_prepares_inputs_and_runs_quota(tmp_path, monkeypatch):
    target_cds = write(tmp_path / "target.fa", ">t1\nATG\n>t2\nATG\n")
    target_gff = write(
        tmp_path / "target.gff",
        "chr1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=t1\nchr1\ttest\tmRNA\t10\t12\t.\t+\t.\tID=t2\n",
    )
    query_cds = write(tmp_path / "query.fa", ">q1\nATG\n>q2\nATG\n")
    query_gff = write(
        tmp_path / "query.gff",
        "chrA\ttest\tmRNA\t1\t3\t.\t+\t.\tID=q1\nchrA\ttest\tmRNA\t10\t12\t.\t+\t.\tID=q2\n",
    )
    commands = []

    def fake_run(command, *, cwd, description):
        assert description
        commands.append(tuple(command))
        write(Path(cwd) / "target.query.lifted.1x2.anchors", "t1\tq1\t10\n")

    monkeypatch.setattr("kffractbias.jcvi._run_checked", fake_run)
    result = run_pairwise_synteny(
        target_cds=target_cds,
        target_gff=target_gff,
        query_cds=query_cds,
        query_gff=query_gff,
        work_dir=tmp_path / "work",
        quota="1:2",
        cpus=2,
        cscore=0.7,
        aligner="last",
    )

    assert result.anchors_path.is_file()
    assert result.target.mapping.matched_gene_count == 2
    assert result.query.mapping.matched_gene_count == 2
    assert "--quota=1:2" in commands[0]
    assert "--cpus=2" in commands[0]


def test_run_pairwise_synteny_rejects_overlapping_ids_before_external_command(
    tmp_path, monkeypatch
):
    target_cds = write(tmp_path / "target.fa", ">shared\nATG\n")
    target_gff = write(tmp_path / "target.gff", "chr1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=shared\n")
    query_cds = write(tmp_path / "query.fa", ">shared\nATG\n")
    query_gff = write(tmp_path / "query.gff", "chrA\ttest\tmRNA\t1\t3\t.\t+\t.\tID=shared\n")

    def unexpected_run(*_args, **_kwargs):
        raise AssertionError("external command must not run")

    monkeypatch.setattr("kffractbias.jcvi._run_checked", unexpected_run)
    with pytest.raises(ValueError, match="must be disjoint"):
        run_pairwise_synteny(
            target_cds=target_cds,
            target_gff=target_gff,
            query_cds=query_cds,
            query_gff=query_gff,
            work_dir=tmp_path / "work",
            quota="1:1",
            cpus=1,
            cscore=0.7,
            aligner="last",
        )


def test_run_self_synteny_uses_native_self_mode_and_symmetric_quota(tmp_path, monkeypatch):
    cds = write(tmp_path / "genome.cds.fa", ">g1\nATG\n>g2\nATG\n")
    gff = write(
        tmp_path / "genome.gff3",
        "##gff-version 3\nchr1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=g1\n"
        "chr2\ttest\tmRNA\t1\t3\t.\t+\t.\tID=g2\n",
    )
    commands = []

    def fake_run(command, *, cwd, description):
        assert description
        commands.append(tuple(command))
        cwd = Path(cwd)
        if command[3] == "align":
            write(cwd / "self.self.last.P98L0.inverse.filtered", "g1\tg2\t95\n")
        else:
            write(cwd / "self.self.lifted.2x2.anchors", "###\ng1\tg2\t10\n")

    monkeypatch.setattr("kffractbias.jcvi._run_checked", fake_run)
    result = run_self_synteny(
        cds=cds,
        gff=gff,
        work_dir=tmp_path / "work",
        depth=2,
        cpus=3,
        cscore=0.8,
        aligner="last",
    )
    assert commands[0][3] == "align"
    assert "--self-hit-percent=98.0" in commands[0]
    assert "--cpus=3" in commands[0]
    assert "kffractbias.selfscan" in commands[1]
    assert "--depth=2" in commands[1]
    assert "--diagonal-bound=300" in commands[1]
    assert result.anchors_path.name == "self.self.lifted.2x2.anchors"
    assert result.depth == 2


def test_run_self_synteny_can_override_jcvi_diagonal_bound(tmp_path, monkeypatch):
    cds = write(tmp_path / "genome.cds.fa", ">g1\nATG\n>g2\nATG\n")
    gff = write(
        tmp_path / "genome.gff3",
        "##gff-version 3\nchr1\ttest\tmRNA\t1\t3\t.\t+\t.\tID=g1\n"
        "chr2\ttest\tmRNA\t1\t3\t.\t+\t.\tID=g2\n",
    )
    commands = []

    def fake_run(command, *, cwd, description):
        assert description
        commands.append(tuple(command))
        cwd = Path(cwd)
        if command[3] == "align":
            write(cwd / "self.self.last.P98L0.inverse.filtered", "g1\tg2\t95\n")
        elif "jcvi.compara.synteny" in command:
            write(cwd / "self.self.anchors", "###\ng1\tg2\t10\n")
            write(cwd / "self.self.lifted.anchors", "###\ng1\tg2\t10\n")
        else:
            write(cwd / "self.self.lifted.1x1.anchors", "###\ng1\tg2\t10\n")

    monkeypatch.setattr("kffractbias.jcvi._run_checked", fake_run)
    result = run_self_synteny(
        cds=cds,
        gff=gff,
        work_dir=tmp_path / "work",
        depth=1,
        cpus=1,
        cscore=0.7,
        aligner="last",
        intrachromosomal_diagonal_bound=7,
    )
    assert len(commands) == 2
    assert "--self-hit-percent=98.0" in commands[0]
    assert "--diagonal-bound=7" in commands[1]
    assert "--depth=1" in commands[1]
    assert result.intrachromosomal_diagonal_bound == 7
