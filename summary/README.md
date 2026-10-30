# Lightweight reference outputs

These files record the manuscript calculations, not new simulation runs.
The inventory CSVs contain 181 sampled times; metadata state the orbit orders.
`inner_runs.json` contains individual scalar and 12-block transport summaries.
`inner_suite.json` groups the two independent seeds per main configuration.
`capture_energy.json` records audited packet-cohort and energy diagnostics.
No full packet records are distributed here. Recreate them with the inner
transport CLI, then use the packet audit to recompute these diagnostics.

The exterior orbit flux is Newtonian time; the inner boundary rest flux includes
the sqrt(0.98) static-to-coordinate clock factor. They should not be equated
without this factor. All summary outputs are specific to the stated preparations.
