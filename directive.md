# Directives

Standing conventions for this campaign, layered on top of `PLAN.md` (the
scientific plan) and `journal.md` (the running log). Check this file for
"how things must be done" rules that apply across many tasks, not just one.
Empty for now — add entries here as decisions get made that should bind
future tasks (mirroring `../Twist/directive.md`'s role for that campaign).

## Engine reuse (superseded 2026-08-22 — see below)

~~Do not fork `../IsingS2/`'s C++ engine (mesh generation, `QfeLatticeS2`,
`QfeIsing`). Build/extend against it in place. If a change there is needed,
make it in `../IsingS2/` and note it in this campaign's `journal.md` with a
pointer to the corresponding `../IsingS2/` commit — do not vendor a copy that
can silently drift.~~

**Reversed 2026-08-22, per explicit user direction**: this package is now
fully self-contained (`src/`, `include/`, `grp/` all vendored copies — see
`CLAUDE.md` "Relationship to `IsingS2/`" and `journal.md` 2026-08-22). The
drift risk the original rule was guarding against is now real and accepted
as a tradeoff, not eliminated: **any future fix landing in `../IsingS2/`
must be manually re-copied here** (no automatic sync) and logged in
`journal.md` with what was copied and why. Before debugging a
mesh/MC-engine issue here, check whether `../IsingS2/` has already fixed it
upstream and whether that fix has been ported into this package's vendored
copy.

## No scratchpad for intermediate/derived data

Never write intermediate or derived files for this campaign (rebinned
jackblocks, temporary analysis outputs, anything downstream of a
production run) to the Claude session scratchpad
(`/scratch/claude-*/.../scratchpad`). It is session-scoped and gets wiped —
anything worth keeping (even a reusable derived artifact, not just final
results) belongs under this repo (e.g. `campaign_runs/<run>/...`) so it
survives past the current conversation and is visible to future sessions.
See `journal.md` 2026-08-23 (push4 jackblocks rebinning) for the incident
that prompted this.

## Closed-form CFT formulas

Never hardcode a from-memory closed-form formula (e.g. the harmonic-level
eigenvalue `C_l(l, Delta)`) into a fitting script without it first being
symbolically derived/verified in a committed script and cross-checked
against known special cases. See `PLAN.md` open items.

## `jack_block_size` must be resized per push, not copied from the last run

Pick `--jack_block_size` so the native block count (`n_traj /
jack_block_size`) lands in the few-hundred range regardless of `n_traj` —
ideally on PLAN.md's own stats gate `n_blocks >= 10*(2*l_max+1)`. Do not
reuse the previous push's literal value: `100` was correct at
`n_traj=20000` (`production_2026-08-23`) but would give 10000 native
blocks/file at `push5_2026-08-24_lmax12`'s `n_traj=1,000,000` — the same
multi-GB-jackblocks-file problem `push4_2026-08-22_lmax12` hit at the
default of 1 (`journal.md` 2026-08-23), just shifted to a different block
size. `push5` used `4000` instead. See `CLAUDE.md` "Cluster conventions".
