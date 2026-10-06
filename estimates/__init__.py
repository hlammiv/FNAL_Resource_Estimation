"""Reproducible resource estimates for the chapters of the FNAL quantum-utility report.

    python -m estimates              # table of every box number, with provenance
    python -m estimates --chapter ch03   # one chapter: assumptions, resources.json rows, exports
    python -m estimates --resources  # regenerate resources.json (feeds the landscape figure)
    python -m estimates --status     # regenerate CIRCUIT_STATUS.md
    python -m pytest estimates/tests # every model vs its chapter's published numbers

One module per chapter, chNN_<name>.py, each exposing:
    Assumptions   frozen dataclass; every input a Tagged value with provenance
    model(a, era) -> Result
    PUBLISHED     {era: Published(...)}  the box numbers as printed, with line refs
See CONTRACT.md.
"""

CHAPTERS = {
    2: "ch02_nu_nucleus",
    3: "ch03_mu2e_0nubb",
    4: "ch04_hybrid_lqcd",
    5: "ch05_qgp_transport",
    6: "ch06_collider",
    7: "ch07_curved_space",
    8: "ch08_baryogenesis",
    9: "ch09_chiral_gauge",
    10: "ch10_finite_density",
}

TEX = {
    2: "applications/app01_neutrino_nucleus.tex",
    3: "applications/app11_mu2e_0nubb.tex",
    4: "applications/app12_hybrid_lqcd.tex",
    5: "applications/app03_qgp_transport.tex",
    6: "applications/app10_collider_physics.tex",
    7: "applications/app04_curved_space.tex",
    8: "applications/app05_baryogenesis.tex",
    9: "applications/app06_chiral_gauge.tex",
    10: "applications/app07_finite_density.tex",
}


def load(ch: int):
    import importlib
    return importlib.import_module(f"estimates.{CHAPTERS[ch]}")


def available() -> list[int]:
    out = []
    for ch in CHAPTERS:
        try:
            load(ch)
            out.append(ch)
        except ImportError:
            pass
    return out
