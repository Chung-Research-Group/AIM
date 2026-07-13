# AIM continuous integration

Pull requests run `.github/workflows/pr-checks.yml`. The workflow is designed to
catch scientific-regression risks without making the existing warning backlog a
merge blocker.

## Required checks

| Check | Purpose | Failure condition |
|---|---|---|
| Repository hygiene | Portable and reviewable repository state | Whitespace errors in changed lines, unresolved conflict markers, case-insensitive path or MATLAB source-name collisions, tracked temporary/build output, or missing CI-critical files |
| MATLAB Code Analyzer (R2024a) | Parse and statically inspect the full `src`, `tests`, and build file | Any Code Analyzer **error**. Warnings are retained in SARIF/MAT artifacts but are not yet blocking |
| MATLAB tests (R2024a / Linux) | Enforce the documented minimum MATLAB release and numerical invariants | Any failed, incomplete, or warning-emitting test |
| MATLAB tests (latest / Windows) | Detect forward-compatibility and primary desktop-platform regressions | Any failed, incomplete, or warning-emitting test |
| PR gate | Stable branch-protection target | Any required job is not successful |

Configure branch protection for `main` to require only the stable check name
**PR gate**. The gate depends on all matrix jobs, so the protected-check list
does not need to change when the matrix expands.

## Local reproduction

From the repository root in MATLAB R2024a or newer:

```matlab
buildtool check
buildtool test
```

The `check` task writes Code Analyzer results to `artifacts/`. The `test` task
writes JUnit, MAT, and Cobertura artifacts to the same directory.

For the non-MATLAB repository checks:

```bash
git diff --check
python ci/check_repository.py
```

## Initial CI baseline

The initial validated workflow run established this baseline:

- 15/15 numerical-kernel tests passed on MATLAB R2024a/Linux and the latest
  MATLAB release/Windows.
- Instrumented line coverage was 197/308 lines (63.96%) across the three
  currently instrumented kernels: `WENO.m`, `Isotherm_functions.m`, and
  `IAST_func_NR.m`.
- Code Analyzer reported no errors. It also reported 100 warnings and 31 notes
  in the current tree; these remain visible in the uploaded SARIF/MAT results.

Coverage is reported but is not yet a merge threshold. A threshold should be
introduced only after end-to-end BreakLab and manuscript golden-case tests are
added, otherwise a percentage target can reward low-value line execution rather
than scientific validation.

## Policy choices

- Code Analyzer errors block immediately; warnings and notes are reported but
  do not yet block. Warning blocking should be introduced after the current
  inventory is triaged and a baseline policy is agreed.
- Numerical tests target identities, limiting cases, symmetry, positivity, and
  finite outputs rather than GUI screenshots.
- Installer compilation is intentionally excluded from required PR checks.
  `AIM_Build.prj` targets a Windows standalone application and MATLAB Compiler
  needs separate CI licensing. Packaging belongs in a manually triggered
  release workflow after a licensing token is configured.
- Third-party GitHub Actions are pinned to the exact commit revisions validated
  by the initial green run. Update those pins in a dedicated dependency PR and
  rerun the full matrix before merging.
- No workflow uses `pull_request_target`; pull-request code executes only with
  read-only repository permissions and without repository secrets.
