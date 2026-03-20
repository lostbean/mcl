# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

mcl is a Haskell implementation of the Markov Cluster Algorithm — a fast, scalable method for unsupervised graph clustering.

## Build & Development

**Prerequisites:** Nix (provides GHC, Cabal, HLS, formatters).

```bash
nix develop                    # Enter dev shell
cabal build --allow-newer      # Build the library
```

## CI & Quality Checks

```bash
bash scripts/ci.sh             # Full CI: format check + build with -Wall -Werror + tests
bash scripts/pre-commit.sh     # Apply formatting to all files
```

CI runs format check, compilation with `-Wall -Werror` (fail on warnings), and tests.

## Formatting

Three formatters via treefmt/nix:
- **Fourmolu** — Haskell (.hs)
- **cabal-fmt** — Cabal files (.cabal)
- **nixpkgs-fmt** — Nix files (.nix)

```bash
nix fmt                    # Format all files
nix fmt -- --ci            # Check formatting without modifying (CI mode)
```

Pre-commit hooks via Lefthook auto-format staged files.

## Git Conventions

- **No footer on commit messages.** Do not add `Co-Authored-By` or any other footer lines.
- Keep commit messages concise and descriptive.
