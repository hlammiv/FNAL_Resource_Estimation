# Releases

A release pairs one commit of the paper repository with one commit of this repository. A pair is valid
when `sync_paper.sh` passes on it: run

```
./sync_paper.sh <path-to-paper-checkout>
```

from a clean checkout of each. The script prints both commit hashes first and marks either tree `-dirty`
if it has uncommitted changes; a dirty run does not certify a pair.

For each posted version of the paper, tag this repository (`vN`, matching the paper's arXiv version) and
add one row below with the two full hashes the script printed.

| version | date | paper commit | models commit |
|---|---|---|---|
