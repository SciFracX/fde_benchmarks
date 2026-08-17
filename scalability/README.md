# Dimension and interval scalability benchmark

This experiment supplements the time-step scaling benchmark in
`memory_fode/`. All three implementations solve the same commensurate system

\[
{}^C D_t^{0.8}y(t)=A_dy(t),
\]

where `A_d` is tridiagonal with `-1.5` on the diagonal and `0.25` on the first
off-diagonals. The initial state is linearly spaced from `1.0` to `0.1`.

The default sweeps isolate two effects:

- Dimension scaling: `d = 5, 10, 15, 20, 25, 30, 35, 40`, with `T = 20`
  and `N = 5120`.
- Interval scaling: `T = 5, 10, 15, 20, 25, 30, 35, 40`, with `d = 10`
  and `h = 1/256`, so `N = 256T`.

Each configuration is run in a fresh process five times. Solver runtime is
measured inside the worker, while peak resident set size (RSS) is measured by
`/usr/bin/time` for the complete Julia, MATLAB, or Python process.

Run the full experiment from the project root:

```bash
julia --project=. scalability/scalability_figure.jl
```

Plot existing results without running workers:

```bash
julia --project=. scalability/scalability_figure.jl --plot-only
```

Results are written to `scalability/scalability_results.csv`; the figure is
saved as PDF, SVG, and PNG in the same directory. The CSV is updated after
every successful worker, so an interrupted experiment can be resumed by
running the same command again.

After both this experiment and `memory_fode/memory_figure.jl` have produced
their CSV files, create the combined time-step/dimension/interval figure with:

```bash
julia --project=. scalability/combined_scalability_figure.jl
```
