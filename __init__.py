"""
hplc_sim.integrations
======================
Thin, dependency-light clients for free/open external resources.
Each module degrades gracefully (returns None / raises a clear error)
if network access or the optional dependency is unavailable, so the
core simulator NEVER hard-depends on any of these.
"""
