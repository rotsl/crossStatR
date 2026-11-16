# Contributing to crossStatR

Thanks for considering a contribution. `crossStatR` is a small package with strong opinions about correctness and reproducibility, so the bar for merging is higher than for a typical utility package. The guidance below explains what the package expects.

## Reporting bugs and requesting features

File issues at <https://github.com/rotsl/crossStatR>. A good bug report contains:

- A minimal reproducible example (a script the package can copy, paste, and run).
- The output of `sessionInfo()`.
- The exact error message, not a paraphrase.
- What the analyst expectsed to happen and what actually happened.

For feature requests, please describe the scientific problem first and the proposed API second. A feature that solves a real analysis problem is more likely to be accepted than one that simply adds a convenience wrapper.

## Pull request workflow

1. Fork the repository and clone the fork.
2. Create a feature branch from `main`. Use a descriptive name: `fix/sdr-perm-seed`, `feat/lod-extra-sigma`, `docs/readme-examples`. Avoid `master`/`main` as branch names for the own work.
3. Make the changes. Keep commits focused; one logical change per commit makes review easier.
4. Run `devtools::check()` locally (see below). The check must pass cleanly (no ERRORs, no WARNINGs, no NOTEs that the analyst cannot justify).
5. Push to the fork and open a pull request against `main`. Describe what the PR does and why, and link any related issues.
6. Respond to review comments. Review timelines depend on maintainer availability.

## Running checks

Before submitting:

```r
# From the package root:
devtools::document()
devtools::test()
devtools::check()
```

`devtools::check()` runs `R CMD check` in a clean temporary library, so it will catch missing `Imports` entries that work locally only because the analyst has the package installed globally. If the change adds a new dependency, add it to `DESCRIPTION` (under `Imports` for runtime, `Suggests` for optional), not just to the function body.

For changes that touch numeric code, also run:

```r
devtools::test_coverage()
```

and aim to keep coverage of new code above 80%.

## Code style

- Use `<-` for assignment, not `=`.
- Indent with two spaces. No tabs.
- Keep lines under 80 characters where reasonable. Long lines make diffs harder to review.
- Name functions in `snake_case`. Internal helpers start with a dot: `.assert_numeric_vector`, `.ci_t`.
- Roxygen comments go above the function definition. Run `devtools::document()` after every change to regenerate `NAMESPACE` and `man/`.
- Match the roxygen conventions in `R/internal_utils.R` and `R/signal_difference_regions.R`: `@param`, `@return`, `@section Methodology:` with `\deqn{}` for display math, `@section Limitations:` (be honest), `@references` (real published works only, do not fabricate), `@examples`, `@export`.

## Tests

The package uses `testthat` (3rd edition). Tests live in `tests/testthat/` with one file per function or per logical group (see `test-catC.R` for the Category C methods).

- Every exported function needs at least one test that exercises a happy path and at least one test that exercises an input-validation failure (an error).
- For numeric code, prefer tests that check properties (e.g. "CI covers the true value in 95% of replicates" using a small simulation) over tests that check exact values, which are fragile.
- Set the seed inside `test_that()` blocks. The global RNG state is not guaranteed.
- Keep simulation sizes small. CI runs the full suite on every push; a 30-second test is fine, a 5-minute test is not.

## Branch and commit conventions

- Branches: `fix/<short-description>`, `feat/<short-description>`, `docs/<short-description>`, `refactor/<short-description>`.
- Commits: imperative mood in the subject line ("Add permutation seed to signal_difference_regions", not "Added" or "Adds"). Subject under 72 characters, body wrapped at 72.

## Scientific content

`crossStatR` documents its methods with explicit math, comparisons to established alternatives, and limitations. If the contribution adds a new method or changes an existing one, the roxygen `@section Methodology:` block must be updated, with display equations where appropriate. Limitations must be honest; do not oversell a method to get it merged.

References in `@references` must be real published works. Do not invent citations to make a method look more established than it is.

## Code of conduct

By participating contributors agree to abide by the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). Reports of misconduct can be sent privately to the maintainer using the contact details in the package metadata. Non-sensitive questions can use the repository issue tracker.
