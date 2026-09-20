#!/usr/bin/env python3
"""Validate SFS, summarize outputs, and plot two method-specific figures.

Only parses outputs already supplied. Does not rerun scientific inference.
Uses Python standard library for parsing; matplotlib is optional for plotting.
"""
import argparse
import csv
import math
from pathlib import Path

CASES = (
    ("Elife", "Elife.YRI.4usfs", "Elife_10samples.final.summary", "Elife.YRI.txt", 20, 5334461, False),
    ("Science", "FitCoal.YRI.usfs", "Science_108samples.final.summary", "Science.YRI.txt", 216, 826650000, False),
    ("MBE", "MBE.YRI.usfs", "MBE_108samples.final.summary", "MBE.YRI.txt", 216, 1350734279, False),
    ("Zhen", "GR.YRI.4fold.sfs", "Zhen_GR_100samples.final.summary", None, 100, 5767108, True),
)
BLOCKBUSTER_MODELS = {
    "Elife": "5_epochs.Elife.yml", "Science": "5_epochs.Science.yml",
    "MBE": "5_epochs.MBE.yml", "Zhen": "4_epochs.Zhen.yml",
}
YEARS_PER_GENERATION = 29
MUTATION_RATE = 1.25e-8


def tsv(path):
    with path.open(newline="") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def sfs_values(path):
    # Uploaded files may contain literal two-character backslash-t, not tabs.
    tokens = path.read_text().replace("\\t", " ").split()
    values = [float(value) for value in tokens]
    if not values or any(not math.isfinite(v) or v < 0 for v in values):
        raise ValueError(f"Invalid nonnegative SFS: {path}")
    return values


def validate_sfs(folder):
    rows = []
    for name, filename, _, _, nseq, length, folded in CASES:
        values = sfs_values(folder / filename)
        if len(values) != nseq - 1:
            raise ValueError(f"{filename}: {len(values)} bins, expected {nseq - 1}")
        if folded and any(values[(nseq // 2):]):
            raise ValueError(f"{filename}: upper folded bins must be zero padded")
        seg = sum(values)
        if seg >= length:
            raise ValueError(f"{filename}: segregating bins exceed total surveyed length")
        rows.append(dict(dataset=name, nseq=nseq, folded=folded, length=length,
                         segregating=round(seg), monomorphic=round(length - seg)))
    return rows


def stairway_minimum(path, lower=500_000, upper=1_500_000):
    # Streaming reads are useful for the ~85,600-row files; duplicate stair steps
    # represent plotting geometry, not independent replicates.
    best = None
    for item in tsv(path):
        year, ne = float(item["year"]), float(item["Ne_median"])
        if lower <= year <= upper and math.isfinite(ne) and ne > 0:
            if best is None or ne < float(best["Ne_median"]):
                best = {k: item[k] for k in ("year", "Ne_median", "Ne_2.5%", "Ne_97.5%")}
    if best is None:
        raise ValueError(f"No Stairway Plot estimates in specified window: {path}")
    return best


def step_intervals(path):
    """Read paired FitCoal breakpoints; use the size on each open interval."""
    rows = [(float(r["year"]), float(r["popSize"])) for r in tsv(path)]
    if any(b[0] < a[0] for a, b in zip(rows, rows[1:])):
        raise ValueError(f"Breakpoints not sorted: {path}")
    return [(a[0], b[0], a[1]) for a, b in zip(rows, rows[1:]) if b[0] > a[0]]


def blockbuster_intervals(path):
    """Parse the best Blockbuster demes model; return (younger, older, Ne)."""
    import yaml
    model = yaml.safe_load(path.read_text())
    if model.get("time_units") != "years" or len(model.get("demes", [])) != 1:
        raise ValueError(f"Expected a single deme with time_units: years: {path}")
    older = float("inf")
    intervals = []
    for epoch in model["demes"][0]["epochs"]:
        younger = float(epoch["end_time"])
        ne = float(epoch["start_size"])
        if not (0 <= younger < older and math.isfinite(ne) and ne > 0):
            raise ValueError(f"Invalid epoch in {path}: {epoch}")
        intervals.append((younger, older, ne))
        older = younger
    if older != 0 or len(intervals) < 2:
        raise ValueError(f"Missing final epoch ending at present: {path}")
    return intervals


def epos_intervals(path):
    """Expand rounded EPOS rows into (dataset, younger year, older year, Ne).

    The reporting convention assigns the size on each row to the interval
    ending at that row's T. Coincident *rounded* endpoints are kept separately
    for annotation, never expanded into an invented duration.
    """
    by_dataset = {}
    for row in tsv(path):
        by_dataset.setdefault(row["dataset"], []).append(row)
    out = {}
    for dataset, rows in by_dataset.items():
        intervals = []
        earlier = 0.0
        for row in rows:
            later = float(row["T_generations"]) * YEARS_PER_GENERATION
            ne = float(row["Ne"])
            if not (later >= earlier and math.isfinite(later)
                    and math.isfinite(ne) and ne > 0):
                raise ValueError(f"Invalid EPOS time or size: {dataset}, {row}")
            intervals.append((earlier, later, ne, row["level"]))
            earlier = later
        out[dataset] = intervals
    return out


def write_tsv(path, records, columns):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", fieldnames=columns)
        writer.writeheader()
        writer.writerows(records)


def summarize(folder, out):
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "sfs_qc.tsv", validate_sfs(folder),
              ["dataset", "nseq", "folded", "length", "segregating", "monomorphic"])
    troughs = []
    for name, _, summary, block, _, _, _ in CASES:
        row = stairway_minimum(folder / summary)
        troughs.append(dict(dataset=name, method="Stairway Plot 2", start_year=row["year"],
                            end_year=row["year"], Ne=row["Ne_median"],
                            Ne_2_5=row["Ne_2.5%"], Ne_97_5=row["Ne_97.5%"],
                            note="minimum median within 0.5-1.5 Ma; bounds are Ne at this time, not a time CI"))
        if block:
            for start, end, ne in step_intervals(folder / block):
                if 500_000 <= start <= 1_500_000 and ne <= 5000:
                    troughs.append(dict(dataset=name, method="FitCoal", start_year=start,
                                        end_year=end, Ne=ne, Ne_2_5="", Ne_97_5="",
                                        note="reported low-size interval (Ne <= 5000); threshold only for display"))
        for young, old, ne in blockbuster_intervals(folder / BLOCKBUSTER_MODELS[name]):
            if young < 1_500_000 and old > 500_000 and ne <= 5000:
                troughs.append(dict(dataset=name, method="Blockbuster", start_year=young,
                                    end_year=old, Ne=ne, Ne_2_5="", Ne_97_5="",
                                    note="best-model epoch; start_year is younger boundary"))
    full_epos = []
    for dataset, intervals in epos_intervals(folder / "epos_reported.tsv").items():
        for young, old, ne, level in intervals:
            note = ("same rounded endpoints; duration unresolved" if young == old
                    else "rounded time boundaries; size assigned to interval ending at T")
            full_epos.append(dict(dataset=dataset, level=level, younger_year=young,
                                  older_year=old, Ne=ne, note=note))
            if level == "3":
                troughs.append(dict(dataset=dataset, method="EPOS",
                                    start_year=young, end_year=old, Ne=ne,
                                    Ne_2_5="", Ne_97_5="", note=note))
    write_tsv(out / "epos_epochs.tsv", full_epos,
              ["dataset", "level", "younger_year", "older_year", "Ne", "note"])
    write_tsv(out / "reported_troughs.tsv", troughs,
              ["dataset", "method", "start_year", "end_year", "Ne", "Ne_2_5", "Ne_97_5", "note"])
    # The mismatch is recorded, not silently corrected; a complete EPOS log is needed.
    (out / "input_caveats.txt").write_text(
        "EPOS MBE log: reported segregating 10,529,304 + monomorphic 1,340,162,791 "
        "= 1,350,692,095, which is 42,184 less than command-line -l 1,350,734,279. "
        "The supplied SFS sums to 10,529,304; cause of mismatch is unknown.\n"
        "FitCoal year/popSize trajectories are supplied for Elife, Science and MBE; "
        "four best-model Blockbuster demes YAML files are supplied, including folded Zhen.\n"
        "The precise provenance of FitCoal-format year/popSize trajectories "
        "should be checked against the original generating logs.\n"
        "EPOS rows were transcribed from the message and are rounded; no machine logs were supplied. "
        "For the full reconstructed curve, each Ne is assigned to the interval ending at its T. "
        "Coincident rounded boundaries are shown as zero-width markers, not invented durations.\n"
    )
    return troughs


def _style(ax, xlim=(.35, 1.6), ylim=(100, 300000)):
    ax.set(xlim=xlim, ylim=ylim, xlabel="Years before present (millions)",
           ylabel="Inferred effective population size (Ne)")
    ax.set_yscale("log")
    ax.grid(alpha=.2)


def plot(folder, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"Elife": "#3B75A6", "Science": "#CC6633", "MBE": "#774499", "Zhen": "#3D8457"}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.8), sharex=True, sharey=True)
    for ax, (name, _, _, fitcoal, _, _, folded) in zip(axes.flat, CASES):
        if fitcoal:
            segments = step_intervals(folder / fitcoal)
            sx = [segments[0][0]/1e6] + [b/1e6 for _, b, _ in segments]
            sy = [n for _, _, n in segments] + [segments[-1][2]]
            ax.step(sx, sy, where="post", color="#284E7D", lw=2, label="FitCoal")
        epochs = blockbuster_intervals(folder / BLOCKBUSTER_MODELS[name])
        # Reverse to ascending time (present -> past) and add vertical jumps.
        bx, by = [], []
        for young, old, ne in reversed(epochs):
            older = min(old, 4_000_000)
            bx.extend([young/1e6, older/1e6])
            by.extend([ne, ne])
        ax.plot(bx, by, color="#C65338", ls="--", lw=2, label="Blockbuster")
        _style(ax, xlim=(.45, 1.3), ylim=(70, 500000))
        ax.set_title(f"{name} YRI" + (" · folded SFS" if folded else ""),
                     loc="left", fontweight="bold")
        if ax not in (axes[0, 0], axes[1, 0]):
            ax.set_ylabel("")
        if ax in (axes[0, 0], axes[0, 1]):
            ax.set_xlabel("")
    axes[1, 1].text(.98, .06, "FitCoal trajectory not supplied",
                    transform=axes[1, 1].transAxes, ha="right", fontsize=8,
                    bbox=dict(facecolor="white", edgecolor="#aaaaaa", alpha=.93, pad=5))
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=2, frameon=False,
               bbox_to_anchor=(.98, 1.005))
    fig.suptitle("Figure 1 | FitCoal and Blockbuster", x=.05, y=.99,
                 ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(out / "figure1_fitcoal_blockbuster.png", dpi=220)
    plt.close(fig)

    # Folded Zhen Stairway Plot intentionally excluded from this figure.
    fig, axes = plt.subplots(1, 3, figsize=(15.3, 5.2), sharex=True, sharey=True)
    epos = epos_intervals(folder / "epos_reported.tsv")
    for ax, (name, _, summary, _, _, _, _) in zip(axes, CASES[:3]):
        xs, ys, last = [], [], None
        for r in tsv(folder / summary):
            point = (float(r["year"]), float(r["Ne_median"]))
            if point[0] > 0 and point != last:
                xs.append(point[0]/1e6)
                ys.append(point[1])
                last = point
        ax.plot(xs, ys, color="#284E7D", lw=1.55,
                label="Stairway Plot 2 median")
        px, py, unresolved = [], [], []
        for young, old, ne, level in epos[name]:
            if old == young:
                unresolved.append((old/1e6, ne))
                continue
            if px and px[-1] == young/1e6:
                px.append(young/1e6)
                py.append(ne)
            px.extend([young/1e6, old/1e6])
            py.extend([ne, ne])
        ax.plot(px, py, color="#B22121", lw=1.9,
                label="EPOS: reconstructed from rounded rows")
        if unresolved:
            ax.scatter(*zip(*unresolved), marker="v", s=90,
                       color="#B22121", edgecolor="white", linewidth=.4,
                       zorder=5, label="EPOS: duration unresolved")
        _style(ax, xlim=(0, 2.65), ylim=(.7, 800000))
        ax.set_title(f"{name} YRI", fontweight="bold")
        if ax is not axes[0]:
            ax.set_ylabel("")
    legend = {}
    for ax in axes:
        for handle, label in zip(*ax.get_legend_handles_labels()):
            legend[label] = handle
    fig.legend(legend.values(), legend.keys(), ncol=3, loc="upper center",
               bbox_to_anchor=(.5, 1.01), frameon=False, fontsize=9)
    fig.suptitle("Figure 2 | Full EPOS trajectories and Stairway Plot 2", x=.01, y=1.06,
                 ha="left", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, .93))
    fig.savefig(out / "figure2_stairway_epos.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data"))
    parser.add_argument("--output", type=Path, default=Path("results"))
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()
    summarize(args.data, args.output)
    if not args.no_plot:
        plot(args.data, args.output)
    print(f"Summaries written to {args.output.resolve()}")


if __name__ == "__main__":
    main()
