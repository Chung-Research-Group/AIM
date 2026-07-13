#!/usr/bin/env python3
"""Fail CI on repository states that are unsafe or non-portable."""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

REQUIRED_PATHS = {
    ".github/workflows/pr-checks.yml",
    "AIM_Build.prj",
    "LICENSE",
    "README.md",
    "buildfile.m",
    "ci/check_repository.py",
    "src/IAST_func_NR.m",
    "src/Isotherm_functions.m",
    "src/Main_app.mlapp",
    "src/WENO.m",
    "tests/test_numerical_kernels.m",
}

TEXT_SUFFIXES = {
    ".gitignore",
    ".gitattributes",
    ".json",
    ".m",
    ".md",
    ".prj",
    ".py",
    ".rst",
    ".txt",
    ".xml",
    ".yaml",
    ".yml",
}

GENERATED_DIRECTORY_NAMES = {
    "__pycache__",
    "aim_build",
    "codegen",
    "slprj",
}

TEMPORARY_FILE_SUFFIXES = (
    ".asv",
    ".autosave",
    ".pyc",
    ".pyo",
    ".swo",
    ".swp",
)

CONFLICT_MARKER = re.compile(r"^(<<<<<<< |=======\s*$|>>>>>>> )")


def tracked_files() -> list[str]:
    result = subprocess.run(
        ["git", "-C", str(ROOT), "ls-files", "-z"],
        check=True,
        capture_output=True,
    )
    return [
        path
        for path in result.stdout.decode("utf-8", errors="surrogateescape").split("\0")
        if path
    ]


def is_text_candidate(path: PurePosixPath) -> bool:
    return path.name.casefold() in TEXT_SUFFIXES or path.suffix.casefold() in TEXT_SUFFIXES


def main() -> int:
    paths = tracked_files()
    path_set = set(paths)
    errors: list[str] = []

    missing = sorted(REQUIRED_PATHS - path_set)
    for path in missing:
        errors.append(f"required path is missing: {path}")

    casefolded_paths: dict[str, list[str]] = defaultdict(list)
    for path in paths:
        casefolded_paths[path.casefold()].append(path)
    for duplicates in casefolded_paths.values():
        if len(duplicates) > 1:
            errors.append(
                "case-insensitive path collision: " + ", ".join(sorted(duplicates))
            )

    source_symbols: dict[str, list[str]] = defaultdict(list)
    for path in paths:
        pure_path = PurePosixPath(path)
        if (
            pure_path.parts
            and pure_path.parts[0].casefold() == "src"
            and pure_path.suffix.casefold() in {".m", ".mlapp"}
        ):
            source_symbols[pure_path.stem.casefold()].append(path)
    for duplicates in source_symbols.values():
        if len(duplicates) > 1:
            errors.append(
                "MATLAB source-name collision on case-insensitive systems: "
                + ", ".join(sorted(duplicates))
            )

    for path in paths:
        pure_path = PurePosixPath(path)
        lower_parts = [part.casefold() for part in pure_path.parts]
        lower_name = pure_path.name.casefold()

        if any(part in GENERATED_DIRECTORY_NAMES for part in lower_parts[:-1]):
            errors.append(f"generated build directory is tracked: {path}")

        if (
            lower_name in {".ds_store", "desktop.ini", "thumbs.db"}
            or lower_name.startswith(".#")
            or lower_name.endswith("~")
            or lower_name.endswith(TEMPORARY_FILE_SUFFIXES)
        ):
            errors.append(f"temporary/editor file is tracked: {path}")

        if not is_text_candidate(pure_path):
            continue

        file_path = ROOT / Path(*pure_path.parts)
        try:
            raw = file_path.read_bytes()
        except OSError as exc:
            errors.append(f"cannot read tracked text file {path}: {exc}")
            continue

        if b"\0" in raw:
            continue

        for line_number, line in enumerate(
            raw.decode("utf-8", errors="replace").splitlines(), start=1
        ):
            if CONFLICT_MARKER.match(line):
                errors.append(
                    f"unresolved merge-conflict marker: {path}:{line_number}"
                )

    if errors:
        print("Repository hygiene check failed:", file=sys.stderr)
        for error in sorted(errors):
            print(f"  - {error}", file=sys.stderr)
        return 1

    print(
        f"Repository hygiene check passed: {len(paths)} tracked paths, "
        f"{len(source_symbols)} MATLAB source names."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
