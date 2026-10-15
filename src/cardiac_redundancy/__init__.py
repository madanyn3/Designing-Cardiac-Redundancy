"""Source hydraulic research model; no clinical or organ-recovery predictions.

Read the package in this order:
1. mechanics: constitutive laws (pressure/flow relationships).
2. source: assemble those laws into the six-reservoir Smith circulation.
3. integration: advance continuous states and switch discrete valve modes.

The runnable entry point is scripts/run_source_baseline.py. A detailed guide is
in docs/m3-code-guide.md. The M2 regional/auxiliary extensions are not yet part
of this package; a successful source run does not establish their feasibility.
"""
