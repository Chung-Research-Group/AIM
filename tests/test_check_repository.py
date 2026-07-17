"""Unit tests for ci/check_repository.py.

The script under test derives its notion of the repository root and the set
of tracked files from module-level state (``ROOT``) and a thin wrapper around
``git ls-files`` (``tracked_files``).  Each test loads a fresh copy of the
module, monkeypatches ``ROOT`` to a temporary directory, and monkeypatches
``tracked_files`` to return a caller-controlled file list. This keeps the
tests hermetic (no real git repository or subprocess calls) while still
exercising the exact validation logic that ships in the script.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "ci" / "check_repository.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "check_repository_under_test", MODULE_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture()
def check_repository(tmp_path, monkeypatch):
    module = load_module()
    monkeypatch.setattr(module, "ROOT", tmp_path)
    return module


def write_tracked_file(root: Path, relative_path: str, content: str = "") -> None:
    file_path = root / relative_path
    file_path.parent.mkdir(parents=True, exist_ok=True)
    file_path.write_text(content, encoding="utf-8")


def populate_minimal_valid_repository(module, root: Path) -> list[str]:
    """Materialize a tracked-file list that satisfies every REQUIRED_PATHS entry.

    Reading the module's own REQUIRED_PATHS (rather than a hand-copied
    literal) keeps this fixture from silently drifting out of sync if the
    script's required-path set changes.
    """
    paths = sorted(module.REQUIRED_PATHS)
    for relative_path in paths:
        write_tracked_file(root, relative_path, content="placeholder content\n")
    return paths


def test_clean_repository_passes(check_repository, tmp_path, monkeypatch, capsys):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 0
    assert "Repository hygiene check passed" in captured.out
    assert captured.err == ""


def test_missing_required_path_is_reported(check_repository, tmp_path, monkeypatch, capsys):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.remove("README.md")
    (tmp_path / "README.md").unlink()
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "required path is missing: README.md" in captured.err


def test_case_insensitive_path_collision_is_reported(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.append("README.MD")
    write_tracked_file(tmp_path, "README.MD", "duplicate\n")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "case-insensitive path collision" in captured.err
    assert "README.md" in captured.err
    assert "README.MD" in captured.err


def test_matlab_source_name_collision_is_reported(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.append("src/weno.m")
    write_tracked_file(tmp_path, "src/weno.m", "% duplicate MATLAB source name\n")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "MATLAB source-name collision" in captured.err
    assert "src/WENO.m" in captured.err
    assert "src/weno.m" in captured.err


def test_matlab_source_and_mlapp_collision_is_reported(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.append("src/Main_app.m")
    write_tracked_file(tmp_path, "src/Main_app.m", "% collides with Main_app.mlapp\n")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "MATLAB source-name collision" in captured.err
    assert "src/Main_app.mlapp" in captured.err
    assert "src/Main_app.m" in captured.err


@pytest.mark.parametrize(
    "generated_path",
    [
        "codegen/report.html",
        "slprj/sim/output.txt",
        "aim_build/AIM.exe",
        "nested/__pycache__/module.pyc",
    ],
)
def test_generated_build_directory_is_reported(
    check_repository, tmp_path, monkeypatch, capsys, generated_path
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.append(generated_path)
    write_tracked_file(tmp_path, generated_path, "generated\n")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert f"generated build directory is tracked: {generated_path}" in captured.err


@pytest.mark.parametrize(
    "temp_path",
    [
        "src/Isotherm_functions.m.asv",
        "notes.txt~",
        ".#lockfile",
        "Thumbs.db",
        "Desktop.ini",
        ".DS_Store",
    ],
)
def test_temporary_or_editor_file_is_reported(
    check_repository, tmp_path, monkeypatch, capsys, temp_path
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.append(temp_path)
    write_tracked_file(tmp_path, temp_path, "temp\n")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert f"temporary/editor file is tracked: {temp_path}" in captured.err


def test_unresolved_conflict_marker_is_reported_with_line_number(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    conflicted_content = (
        "line one\n"
        "<<<<<<< HEAD\n"
        "our change\n"
        "=======\n"
        "their change\n"
        ">>>>>>> feature-branch\n"
    )
    write_tracked_file(tmp_path, "README.md", conflicted_content)
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "unresolved merge-conflict marker: README.md:2" in captured.err
    assert "unresolved merge-conflict marker: README.md:6" in captured.err
    # The bare '=======' separator must not be misclassified as a marker.
    assert "README.md:4" not in captured.err


def test_markdown_heading_of_equals_signs_is_not_a_false_positive(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    write_tracked_file(tmp_path, "README.md", "Title\n========\nBody text\n")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    assert exit_code == 0
    assert capsys.readouterr().err == ""


def test_binary_file_with_null_byte_is_skipped_for_conflict_scan(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    binary_path = tmp_path / "README.md"
    binary_path.write_bytes(b"binary\x00<<<<<<< HEAD\x00garbage")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    # Presence of a NUL byte causes the conflict-marker scan to be skipped
    # entirely for that file, rather than raising a decode error.
    assert exit_code == 0
    assert capsys.readouterr().err == ""


def test_unreadable_tracked_text_file_is_reported(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.append("MISSING.md")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    assert "cannot read tracked text file MISSING.md" in captured.err


def test_multiple_errors_are_all_reported_and_sorted(
    check_repository, tmp_path, monkeypatch, capsys
):
    paths = populate_minimal_valid_repository(check_repository, tmp_path)
    paths.remove("LICENSE")
    (tmp_path / "LICENSE").unlink()
    paths.append("Thumbs.db")
    write_tracked_file(tmp_path, "Thumbs.db", "temp\n")
    monkeypatch.setattr(check_repository, "tracked_files", lambda: paths)

    exit_code = check_repository.main()

    captured = capsys.readouterr()
    assert exit_code == 1
    lines = [line for line in captured.err.splitlines() if line.strip().startswith("-")]
    assert lines == sorted(lines)
    assert any("required path is missing: LICENSE" in line for line in lines)
    assert any("temporary/editor file is tracked: Thumbs.db" in line for line in lines)


@pytest.mark.parametrize(
    "name,expected",
    [
        ("README.md", True),
        ("script.py", True),
        (".gitignore", True),
        (".gitattributes", True),
        ("data.bin", False),
        ("image.png", False),
        ("archive.tar.gz", False),
    ],
)
def test_is_text_candidate(check_repository, name, expected):
    from pathlib import PurePosixPath

    assert check_repository.is_text_candidate(PurePosixPath(name)) is expected