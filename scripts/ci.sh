#!/usr/bin/env bash
set -euo pipefail

echo "=== Format Check ==="
nix fmt -- --ci

echo ""
echo "=== Build (-Wall -Werror) ==="
nix develop --command cabal build --allow-newer --ghc-options="-Wall -Werror"

echo ""
echo "=== Tests (-Wall -Werror) ==="
nix develop --command cabal test --allow-newer --ghc-options="-Wall -Werror"

echo ""
echo "=== CI Passed ==="
