#!/bin/sh
# KSHELL interaction files (KB3G: Poves et al., NPA 694, 157 (2001); GXPF1A: Honma et al., EPJA 25 s01, 499 (2005);
# USDB: Brown and Richter, PRC 74, 034315 (2006); the KSHELL file header says USDA but the single-particle
# energies are USDB's), as distributed with KSHELL (github.com/GaffaSnobb/kshell), pinned to one commit and
# checked by sha256. Not mirrored here; fetch before running kb3g_lam.py, run_gaps.py or al27_sd.py.
set -e
COMMIT=41ea03a6707fb588e176e906edea85be2a46ac79
DIR="$(dirname "$0")/snt"
mkdir -p "$DIR"
for f in kb3g gxpf1a usdb; do
  curl -sf -o "$DIR/$f.snt" "https://raw.githubusercontent.com/GaffaSnobb/kshell/$COMMIT/snt/$f.snt"
done
cat > "$DIR/SHA256SUMS" <<EOF
ca5efb317f3ae5ec259337dec5a7132c269c1064ce9f47e568a5bb9ccc0e1142  kb3g.snt
3feac25ea63a087925a63f5b0afe7b6c54fbac67e21ff5b3dd16519fc06f1c9a  gxpf1a.snt
4e70131746d4aa2d90eb352f99c32744a140468a3b5c36663fa825a139f866e6  usdb.snt
EOF
(cd "$DIR" && sha256sum -c SHA256SUMS)
