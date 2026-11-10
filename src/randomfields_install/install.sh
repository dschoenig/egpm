#!/usr/bin/env bash
# Build RandomFieldsUtils 1.2.5 and RandomFields 3.3.14 from the CRAN archive
# on current R. Tested with R 4.6.1 and gcc 13.3 on Ubuntu 24.04.
#
# Usage: bash install.sh [library_dir]
# Without library_dir, installs into the first entry of .libPaths().
#
# On an Alliance cluster, run it on a login node (compute nodes have no
# internet), after loading the same modules the jobs use, e.g.
#   module load StdEnv/2026 gcc/15.3.0 r/4.6.1
# The installed packages only work with that R minor version (4.6.x).
#
# The patch stops a segfault in every RFsimulate() call: option parsing read
# the first element of an empty vector, and R returns the invalid pointer 0x1
# for the data of an empty vector.

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if ! command -v R >/dev/null; then
  echo "R not found; load the R module first (module load r/...)." >&2
  exit 1
fi

lib="${1:-$(Rscript -e 'cat(.libPaths()[1])')}"
mkdir -p "$lib"
lib="$(cd "$lib" && pwd)"
echo "Installing into $lib with $(R --version | head -n 1)"

archive="https://cran.r-project.org/src/contrib/Archive"
build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

fetch() {
  if ! curl -fsS -m 120 -o "$2" "$1"; then
    echo "Download failed: $1" >&2
    echo "Compute nodes usually have no internet; run this on a login node." >&2
    exit 1
  fi
}
fetch "$archive/RandomFieldsUtils/RandomFieldsUtils_1.2.5.tar.gz" "$build/rfu.tar.gz"
fetch "$archive/RandomFields/RandomFields_3.3.14.tar.gz" "$build/rf.tar.gz"

tar -xzf "$build/rfu.tar.gz" -C "$build"
patch -p1 -d "$build" -i "$here/RandomFieldsUtils_1.2.5_empty-vector.patch"

export R_MAKEVARS_USER="$here/Makevars"
export R_LIBS="$lib${R_LIBS:+:$R_LIBS}"
export MAKEFLAGS="${MAKEFLAGS:--j4}"

R CMD INSTALL --no-docs --library="$lib" "$build/RandomFieldsUtils"
R CMD INSTALL --no-docs --library="$lib" "$build/rf.tar.gz"

Rscript -e "library(RandomFields, lib.loc = '$lib'); x <- RFsimulate(RMexp(), x = 1:10); cat('RFsimulate OK\n')"
