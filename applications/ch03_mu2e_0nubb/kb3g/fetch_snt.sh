#!/bin/sh
# KSHELL interaction files (KB3G: Poves et al., NPA 694, 157 (2001); GXPF1A: Honma et al., EPJA 25 s01, 499 (2005)),
# as distributed with KSHELL (github.com/GaffaSnobb/kshell). Not mirrored here; fetch before running kb3g_lam.py.
mkdir -p "$(dirname "$0")/snt"
for f in kb3g gxpf1a; do
  curl -s -o "$(dirname "$0")/snt/$f.snt" "https://raw.githubusercontent.com/GaffaSnobb/kshell/master/snt/$f.snt"
done
ls -l "$(dirname "$0")/snt"
