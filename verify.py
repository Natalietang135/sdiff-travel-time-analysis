"""Check saved results against separate ray-path calls and missing-arrival cases."""

import csv
from collections import Counter
import math

from run import OUT, DEPTHS_KM, DISTANCES_DEG, PHASES, TauPyModel

rows = list(csv.DictReader((OUT / "arrivals.csv").open()))
model = TauPyModel("prem")
expected = {(str(h), str(d), p) for h in DEPTHS_KM for d in DISTANCES_DEG for p in PHASES}
assert {(r["source_depth_km"], r["distance_deg"], r["phase"]) for r in rows} == expected
for r in rows:
    if r["time_s"]:
        assert math.isfinite(float(r["time_s"])) and float(r["time_s"]) > 0
    if r["status"] == "ok":
        delta = float(r["time_s"]) - float(r["sdiff_time_s"])
        assert abs(delta - float(r["relative_time_s"])) < 1e-9
        assert r["in_0_40_s"] == str(0 <= delta <= 40)
    else:
        assert r["relative_time_s"] == r["in_0_40_s"] == ""
        if r["status"] == "no_phase_arrival":
            assert r["time_s"] == r["arrival_index"] == ""
        else:
            assert r["status"] == "no_sdiff_reference" and r["sdiff_time_s"] == ""

max_error = 0
checked = 0
for depth, distance in [(20, 90), (20, 105), (100, 110), (300, 110), (300, 120), (600, 130), (600, 150)]:
    selected = [r for r in rows if int(r["source_depth_km"]) == depth and
                int(r["distance_deg"]) == distance]
    for phase in PHASES:
        actual = [r for r in selected if r["phase"] == phase and r["time_s"]]
        rays = model.get_ray_paths(depth, distance, phase_list=[phase])
        assert len(actual) == len(rays), (depth, distance, phase)
        for row, ray in zip(actual, rays):
            error = abs(float(row["time_s"]) - ray.time)
            max_error = max(max_error, error)
            assert error < 0.01, (depth, distance, phase, error)
            assert abs(float(row["purist_distance_deg"]) - ray.purist_distance) < 1e-6
            assert abs(ray.path[-1]["time"] - ray.time) < 0.01
            checked += 1

branches = [r for r in rows if r["source_depth_km"] == "300" and
            r["distance_deg"] == "110" and r["phase"] == "SKKS"]
assert [round(float(r["purist_distance_deg"])) for r in branches] == [110, 250]
messages = [
    f"PASS: {len(expected)} depth/distance/phase combinations represented; {len(rows)} rows.",
    f"PASS: all relative times, inclusive window flags and missing values consistent.",
    f"PASS: {checked} arrivals checked with separate single-phase get_ray_paths calls.",
    f"Maximum travel-time difference: {max_error:.8f} s (tolerance 0.01 s).",
    "PASS: SKKS 110-degree and 250-degree paths both retained at 300 km / 110 degrees.",
    f"Status counts: {dict(Counter(r['status'] for r in rows))}",
    "Scope: internal numerical checks using ObsPy; not an independent Earth-model validation.",
]
(OUT / "verification.txt").write_text("\n".join(messages) + "\n")
print("\n".join(messages))
