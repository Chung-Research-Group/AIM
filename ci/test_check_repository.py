#!/usr/bin/env python3
"""Unit tests for ci/check_repository.py.

These tests exercise the repository-hygiene checker in isolation by
pointing its ``ROOT`` module attribute at a disposable temporary
directory and stubbing ``tracked_files`` so real Git state is never
required. Only ``main`` (and the small pure helpers it calls) is
exercised here; no test depends on the actual contents of this
repository, except for the end-to-end regression test at the bottom of
the file which intentionally runs the script against the real
repository to confirm it still reports success.
"""

from __future__ import annotations

import importlib.util
import io
import subprocess
import sys
import tempfile
import unittest
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from pathlib import Path
from unittest import mock

MODULE_PATH = Path(__file__).resolve().with_name("check_repository.py")
_SPEC = importlib.util.spec_from_file_location(
    "check_repository_under_test", MODULE_PATH
)
check_repository = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(check_repository)


def _write_repo(root: Path, paths: set[str], content: bytes = b"placeholder\n") -> None:
    """Materialize every path in ``paths`` under ``root`` with ``content``."""
    for relative_path in paths:
        file_path = root / Path(relative_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_bytes(content)


@contextmanager
def _temp_repo_root():
    with tempfile.TemporaryDirectory() as tmp_dir:
        yield Path(tmp_dir)


def _run_main(root: Path, paths: set[str]) -> tuple[int, str, str]:
    """Run check_repository.main() against a fake ROOT/tracked-file set."""
    stdout = io.StringIO()
    stderr = io.StringIO()
    with mock.patch.object(check_repository, "ROOT", root), mock.patch.object(
        check_repository, "tracked_files", return_value=sorted(paths)
    ), redirect_stdout(stdout), redirect_stderr(stderr):
        return_code = check_repository.main()
    return return_code, stdout.getvalue(), stderr.getvalue()


class RequiredPathsTests(unittest.TestCase):
    def test_couple_pressure_temperature_rates_module_is_required(self):
        # Regression check for this PR: the new coupling module must be
        # listed so CI fails if it is ever accidentally deleted.
        self.assertIn(
            "src/couple_pressure_temperature_rates.m",
            check_repository.REQUIRED_PATHS,
        )


class IsTextCandidateTests(unittest.TestCase):
    def test_matches_known_suffix(self):
        self.assertTrue(
            check_repository.is_text_candidate(check_repository.PurePosixPath("src/foo.py"))
        )

    def test_matches_dotfile_by_full_name(self):
        self.assertTrue(
            check_repository.is_text_candidate(
                check_repository.PurePosixPath(".gitignore")
            )
        )

    def test_is_case_insensitive(self):
        self.assertTrue(
            check_repository.is_text_candidate(check_repository.PurePosixPath("FOO.PY"))
        )

    def test_rejects_unknown_suffix(self):
        self.assertFalse(
            check_repository.is_text_candidate(check_repository.PurePosixPath("image.png"))
        )

    def test_rejects_extensionless_license_file(self):
        # LICENSE has no suffix and its lowercase name is not in
        # TEXT_SUFFIXES, so it is intentionally never content-scanned.
        self.assertFalse(
            check_repository.is_text_candidate(check_repository.PurePosixPath("LICENSE"))
        )

    def test_rejects_mlapp_binary_container(self):
        self.assertFalse(
            check_repository.is_text_candidate(
                check_repository.PurePosixPath("src/Main_app.mlapp")
            )
        )


class TrackedFilesTests(unittest.TestCase):
    def test_splits_nul_delimited_output_and_drops_empty_entries(self):
        completed = subprocess.CompletedProcess(
            args=["git"], returncode=0, stdout=b"a.py\0dir/b.py\0"
        )
        with mock.patch.object(
            check_repository.subprocess, "run", return_value=completed
        ) as run_mock:
            result = check_repository.tracked_files()

        self.assertEqual(result, ["a.py", "dir/b.py"])
        run_mock.assert_called_once()
        called_args = run_mock.call_args.args[0]
        self.assertEqual(called_args[:2], ["git", "-C"])
        self.assertIn("ls-files", called_args)
        self.assertIn("-z", called_args)


class MainRepositoryChecksTests(unittest.TestCase):
    def _valid_paths(self) -> set[str]:
        return set(check_repository.REQUIRED_PATHS)

    def test_passes_on_a_minimal_valid_repository(self):
        with _temp_repo_root() as root:
            paths = self._valid_paths()
            _write_repo(root, paths)

            return_code, stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 0, stderr)
        self.assertIn("Repository hygiene check passed", stdout)
        self.assertEqual(stderr, "")

    def test_reports_missing_required_path(self):
        with _temp_repo_root() as root:
            paths = self._valid_paths()
            paths.discard("src/couple_pressure_temperature_rates.m")
            _write_repo(root, paths)

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn(
            "required path is missing: src/couple_pressure_temperature_rates.m",
            stderr,
        )

    def test_detects_case_insensitive_path_collision(self):
        with _temp_repo_root() as root:
            paths = {"docs/Notes.txt", "docs/notes.txt"}
            _write_repo(root, paths)

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn("case-insensitive path collision", stderr)
        self.assertIn("docs/Notes.txt", stderr)
        self.assertIn("docs/notes.txt", stderr)

    def test_detects_matlab_source_name_collision(self):
        with _temp_repo_root() as root:
            paths = {"src/Foo.m", "src/foo.mlapp"}
            _write_repo(root, paths)

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn(
            "MATLAB source-name collision on case-insensitive systems", stderr
        )

    def test_ignores_matlab_source_collision_outside_src(self):
        with _temp_repo_root() as root:
            paths = {"tests/Foo.m", "examples/foo.m"}
            _write_repo(root, paths)

            _return_code, _stdout, stderr = _run_main(root, paths)

        self.assertNotIn("MATLAB source-name collision", stderr)

    def test_detects_generated_directory_tracked(self):
        with _temp_repo_root() as root:
            paths = {"codegen/generated_kernel.m"}
            _write_repo(root, paths)

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn(
            "generated build directory is tracked: codegen/generated_kernel.m",
            stderr,
        )

    def test_detects_temporary_editor_file_tracked(self):
        with _temp_repo_root() as root:
            paths = {"module.pyc"}
            # No file is created on disk: the temporary-file check is
            # purely name-based and must not require reading the file.

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn("temporary/editor file is tracked: module.pyc", stderr)

    def test_detects_ds_store_tracked(self):
        with _temp_repo_root() as root:
            paths = {"src/.DS_Store"}

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn("temporary/editor file is tracked", stderr)

    def test_detects_unresolved_merge_conflict_marker(self):
        with _temp_repo_root() as root:
            paths = {"docs/notes.md"}
            _write_repo(
                root,
                paths,
                content=b"before\n<<<<<<< HEAD\nmine\n=======\ntheirs\n>>>>>>> branch\nafter\n",
            )

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn("unresolved merge-conflict marker: docs/notes.md:2", stderr)
        self.assertIn("unresolved merge-conflict marker: docs/notes.md:6", stderr)

    def test_equals_sign_only_line_is_not_a_conflict_marker(self):
        # A standalone '=======' line is valid Markdown/RST and must not
        # be misreported as an unresolved conflict marker.
        with _temp_repo_root() as root:
            paths = self._valid_paths() | {"docs/notes.md"}
            _write_repo(root, self._valid_paths())
            (root / "docs" / "notes.md").parent.mkdir(parents=True, exist_ok=True)
            (root / "docs" / "notes.md").write_bytes(b"Title\n=======\nBody text.\n")

            return_code, stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 0, stderr)
        self.assertIn("Repository hygiene check passed", stdout)

    def test_skips_content_scan_for_files_containing_nul_bytes(self):
        # Binary files that happen to carry a text-like suffix must be
        # skipped entirely, even if their bytes contain a conflict-marker
        # look-alike sequence.
        with _temp_repo_root() as root:
            paths = self._valid_paths() | {"src/binary_blob.m"}
            _write_repo(root, self._valid_paths())
            (root / "src" / "binary_blob.m").write_bytes(
                b"\x00binary payload\n<<<<<<< HEAD\n"
            )

            return_code, stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 0, stderr)
        self.assertIn("Repository hygiene check passed", stdout)

    def test_reports_unreadable_tracked_text_file(self):
        with _temp_repo_root() as root:
            paths = {"src/missing_on_disk.m"}
            # Intentionally do not create the file.

            return_code, _stdout, stderr = _run_main(root, paths)

        self.assertEqual(return_code, 1)
        self.assertIn("cannot read tracked text file src/missing_on_disk.m", stderr)


class RealRepositoryEndToEndTest(unittest.TestCase):
    """Runs the unmodified script against the actual checked-out repository.

    This complements the isolated unit tests above by confirming the
    real script (invoked exactly as CI does) still succeeds now that
    ``src/couple_pressure_temperature_rates.m`` has been added.
    """

    def test_real_repository_passes_hygiene_check(self):
        result = subprocess.run(
            [sys.executable, str(MODULE_PATH)],
            cwd=MODULE_PATH.parents[1],
            capture_output=True,
            text=True,
            check=False,
        )

        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("Repository hygiene check passed", result.stdout)


if __name__ == "__main__":
    unittest.main()