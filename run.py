"""PREM timing exercise. Run with python run.py; see README.md."""

import csv
import json
import os
from pathlib import Path
import platform
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sdiff-matplotlib"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import obspy
from obspy.taup import TauPyModel

# Exploratory parameter choices. Receiver is at the surface.
MODEL = "prem"
DEPTHS_KM = [20, 100, 300, 600]
DISTANCES_DEG = range(90, 151)
PHASES = ["Sdiff", "SKS", "SKKS", "pSKS", "sSKS", "pSKKS", "sSKKS"]
WINDOW_S = (0, 40)
ZOOM_S = (-120, 180)  # Display limits only; no pre-arrival screening boundary.
OUT = Path(__file__).resolve().parent / "output"
COLORS = ["#111111", "#0072B2", "#D55E00", "#009E73", "#56B4E9", "#CC79A7", "#8C6500"]
MARKERS = [".", "o", "s", "^", "v", "D", "P"]


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def grid_ranges(values):
    """Compress consecutive integer grid samples; these are not exact crossings."""
    groups = []
    for value in sorted(set(values)):
        if groups and value == groups[-1][-1] + 1:
            groups[-1].append(value)
        else:
            groups.append([value])
    return "; ".join(str(g[0]) if len(g) == 1 else f"{g[0]}-{g[-1]}" for g in groups)


def calculate(model):
    rows, summary = [], []
    for depth in DEPTHS_KM:
        by_distance = {}
        for distance in DISTANCES_DEG:
            arrivals = model.get_travel_times(
                source_depth_in_km=depth, distance_in_degree=distance,
                phase_list=PHASES, receiver_depth_in_km=0.0,
            )
            by_phase = {phase: [a for a in arrivals if a.name == phase] for phase in PHASES}
            by_distance[distance] = by_phase
            refs = by_phase["Sdiff"]
            for phase in PHASES:
                # Keep every arrival and every Sdiff reference, including major arcs.
                for index, arrival in enumerate(by_phase[phase] or [None], 1):
                    for ref_index, ref in enumerate(refs or [None], 1):
                        relative = float(arrival.time - ref.time) if arrival and ref else None
                        rows.append(dict(
                            model=MODEL, source_depth_km=depth, distance_deg=distance,
                            phase=phase, arrival_index=index if arrival else None,
                            time_s=float(arrival.time) if arrival else None,
                            purist_distance_deg=float(arrival.purist_distance) if arrival else None,
                            ray_param_s_per_deg=float(arrival.ray_param_sec_degree) if arrival else None,
                            sdiff_reference_index=ref_index if ref else None,
                            sdiff_time_s=float(ref.time) if ref else None,
                            relative_time_s=relative,
                            status="no_phase_arrival" if arrival is None else
                                   "no_sdiff_reference" if ref is None else "ok",
                            in_0_40_s=(WINDOW_S[0] <= relative <= WINDOW_S[1])
                                      if relative is not None else None,
                        ))
        for phase in PHASES:
            valid = [r for r in rows if r["source_depth_km"] == depth and
                     r["phase"] == phase and r["relative_time_s"] is not None]
            closest = min(valid, key=lambda r: abs(r["relative_time_s"] - 40), default=None)
            summary.append(dict(
                model=MODEL, source_depth_km=depth, phase=phase,
                no_arrival_grid_points=sum(not p[phase] for p in by_distance.values()),
                multiple_arrival_grid_points=sum(len(p[phase]) > 1 for p in by_distance.values()),
                sdiff_available_distance_samples_deg=grid_ranges(
                    [d for d, p in by_distance.items() if p["Sdiff"]]),
                in_0_40_distance_samples_deg=grid_ranges(
                    [r["distance_deg"] for r in valid if r["in_0_40_s"]]),
                closest_to_40_distance_deg=closest["distance_deg"] if closest else None,
                closest_to_40_relative_time_s=closest["relative_time_s"] if closest else None,
            ))
        print(f"Calculated depth {depth} km ({len(DISTANCES_DEG)} distances)", flush=True)
    return rows, summary


def plot(rows, relative):
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), sharex=True, sharey=True)
    key = "relative_time_s" if relative else "time_s"
    for ax, depth in zip(axes.flat, DEPTHS_KM):
        selected = [r for r in rows if r["source_depth_km"] == depth]
        if relative:
            ax.axhspan(*WINDOW_S, color="#F3D46B", alpha=0.28)
            ax.axhline(0, color="0.4", lw=0.7)
            ax.axhline(40, color="#947000", lw=0.8, ls="--")
            ax.set_ylim(*ZOOM_S)
            absent = sorted({r["distance_deg"] for r in selected if r["sdiff_time_s"] is None})
            ax.scatter(absent, [0.025] * len(absent), transform=ax.get_xaxis_transform(),
                       marker="x", color="0.6", s=14, clip_on=False)
        for phase, color, marker in zip(PHASES, COLORS, MARKERS):
            points = [r for r in selected if r["phase"] == phase and r[key] is not None]
            # Scatter avoids inventing connections between branches or across gaps.
            ax.scatter([r["distance_deg"] for r in points], [r[key] for r in points],
                       color=color, marker=marker, s=12, linewidths=0, label=phase)
        ax.set_title(f"Source depth {depth} km", loc="left", fontweight="bold")
        ax.set_xlim(min(DISTANCES_DEG) - 1, max(DISTANCES_DEG) + 1)
        ax.grid(alpha=0.18)
        ax.set_axisbelow(True)
        ax.set_xlabel("Epicentral distance (degrees)")
        ax.set_ylabel("Time relative to Sdiff (s)" if relative else "Travel time from origin (s)")
    fig.suptitle("PREM | Phases relative to Sdiff" if relative else
                 "PREM | Travel times, including all returned arrivals", fontsize=17, y=0.985)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=7, frameon=False,
               bbox_to_anchor=(0.5, 0.947), markerscale=2)
    caption = ("Shading: provisional 0-40 s screen. Grey x: no computed Sdiff reference.\n"
               "Vertical limits are a display zoom only; all arrivals remain in arrivals.csv and the travel-time plot."
               if relative else
               "20, 100, 300, 600 km; 90-150 degrees in 1-degree steps; surface receiver.\n"
               "Each point is a computed arrival. Multiple branches and major-arc paths are retained; gaps are not zero travel times.")
    fig.text(0.5, 0.025, caption, ha="center", fontsize=10, color="0.3")
    fig.tight_layout(rect=(0, 0.08, 1, 0.89))
    fig.savefig(OUT / ("relative_times.png" if relative else "travel_times.png"), dpi=180)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    model = TauPyModel(model=MODEL)
    rows, summary = calculate(model)
    write_csv("arrivals.csv", rows)
    write_csv("screening_summary.csv", summary)
    metadata = dict(
        earth_model=MODEL, python_version=platform.python_version(),
        obspy_version=obspy.__version__, matplotlib_version=matplotlib.__version__,
        source_depths_km=DEPTHS_KM, distances_deg=list(DISTANCES_DEG),
        receiver_depth_km=0, phases=PHASES, provisional_window_s=WINDOW_S,
        relative_plot_limits_s=ZOOM_S, ray_param_tol_s=0.1,
        assumptions="Exploratory parameter choices; see README.md.",
        relative_reference="Every returned Sdiff arrival at the same depth/distance; all pairs retained.",
        arrival_index="Time order within one phase/query, not a persistent branch identifier.",
        missing_values="Blank numeric cells are missing, never zero. Check status.",
        screening="Inclusive 0 <= phase time - Sdiff time <= 40 s; not an amplitude or safety prediction.",
        plot_style="Unconnected points preserve gaps and avoid invented branch connections.",
    )
    (OUT / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    plot(rows, relative=False)
    plot(rows, relative=True)
    print(f"Saved {len(rows)} arrival/reference rows and figures to {OUT}")


if __name__ == "__main__":
    main()
