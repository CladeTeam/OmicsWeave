"""Draw the OmicsWeave logo: assets/logo.svg and assets/logo-animated.svg.

The mark is a single thread that circles twice before it closes, so its two
laps read as the complementary strands of a double helix bent into a ring
(the torus knot T(2, 2K) for half-integer K). Its colour cycles through the
four base colours once over both laps, so the two laps at any angle sit half
a cycle apart: A faces T and C faces G.

Over/under is drawn by removing the lower pass wherever it comes within GAP
of the upper pass. Nothing is painted in a background colour, so the rungs
stay visible and the logo has a transparent background.

Colour depends only on polar angle, so it is painted as a static fan of
wedges and each thread segment is cut out of its fan by a mask. The animated
logo turns the masks and rungs while the fans stay put: the crossings travel
round the ring and the colours hold still, as if the helix were spinning
about its own axis. Turning by two crossings maps the geometry onto itself,
so that spin loops seamlessly; the whole mark also turns slowly.

Usage: python scripts/plot_logo.py [--out-dir DIR]
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

# Thread: radius r(t) = R + AMP * sin(K * t) at polar angle t, t in [0, 4*pi).
K = 3.5  # half-integer: one closed thread with 2K (odd) crossings
R = 80.0
AMP = 18.0
STROKE = 13.0
GAP = 4.5  # clearance between the upper pass and the cut ends of the lower one
MARGIN = 2.0

# Colours, cycled along the thread. Half a cycle apart: A<->T, C<->G.
A, C, G, T = "#1F9E89", "#3B5BA9", "#E8A33D", "#E0604D"
CYCLE = (A, C, T, G)
RUNG = "#C9D2DC"

# Rungs are radial sectors: rung and gap widen together towards the rim, so
# their ratio is the same at every radius.
RUNGS_PER_LOBE = 5
RUNG_PITCH = 0.08  # radians between rung centres
RUNG_DUTY = 0.5  # fraction of the pitch covered by a rung
RUNG_TUCK = STROKE / 4  # rung ends stop this far inside the strand centre line
RUNG_MIN_VISIBLE = 3.0  # drop rungs whose exposed length is shorter than this

# Colour fans: overlapping wedges; the masks cut the only visible edges.
WEDGE = math.radians(1.5)
WEDGE_OVERLAP = 0.5  # fraction of a wedge each one extends into the next
FAN_PAD = 0.1  # radians of fan beyond each end of a segment
CUT_REACH = 0.45  # cut length either side of a crossing, as a fraction of DT

# Animation: seconds for the helix to advance by two crossings, and for the
# whole mark to make one turn (0 keeps it upright). Positive is clockwise.
SPIN_PERIOD = 4.0
TURN_PERIOD = 40.0

CROSSINGS = round(2 * K)
PASSES = 2 * CROSSINGS  # each crossing is visited twice along the thread
DT = math.pi / K  # parameter step between consecutive crossing passes
SPIN = 2 * DT  # smallest rotation that maps segment i onto segment i + 1
SIZE = 2 * (R + AMP + STROKE / 2 + MARGIN)
CENTER = SIZE / 2


def fmt(v: float) -> str:
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if s == "-0" else s


def point_at(t: float, r: float) -> tuple[float, float]:
    """Point at polar angle t (clockwise from 12 o'clock) and radius r."""
    return CENTER + r * math.sin(t), CENTER - r * math.cos(t)


def point(t: float) -> tuple[float, float]:
    return point_at(t, R + AMP * math.sin(K * t))


def polyline(pts: list[tuple[float, float]]) -> str:
    return "M" + " L".join(f"{fmt(x)},{fmt(y)}" for x, y in pts)


def thread(t0: float, t1: float, steps: int) -> str:
    return polyline([point(t0 + (t1 - t0) * i / steps) for i in range(steps + 1)])


def colour(t: float) -> str:
    """Thread colour at parameter t: one full cycle over both laps."""
    x = (t / (4 * math.pi)) % 1.0 * len(CYCLE)
    i = int(x)
    f = x - i
    c1, c2 = CYCLE[i], CYCLE[(i + 1) % len(CYCLE)]
    return "#" + "".join(
        f"{round(int(c1[j : j + 2], 16) * (1 - f) + int(c2[j : j + 2], 16) * f):02X}"
        for j in (1, 3, 5)
    )


def rotation(period: float, degrees: float) -> str:
    """SMIL rotation about the centre, looping forever; empty when disabled."""
    if period <= 0:
        return ""
    c = fmt(CENTER)
    return (
        f'<animateTransform attributeName="transform" type="rotate" '
        f'from="0 {c} {c}" to="{fmt(degrees)} {c} {c}" dur="{fmt(period)}s" '
        f'repeatCount="indefinite"/>'
    )


def rung(t0: float, t1: float, steps: int = 4) -> str:
    """Sector between angles t0 and t1, spanning the two strands.

    Its ends follow the strand centre lines, so they stay hidden under the
    strands drawn on top."""
    ts = [t0 + (t1 - t0) * i / steps for i in range(steps + 1)]
    half_span = [AMP * abs(math.sin(K * t)) - RUNG_TUCK for t in ts]
    inner = [point_at(t, R - s) for t, s in zip(ts, half_span, strict=True)]
    outer = [point_at(t, R + s) for t, s in zip(reversed(ts), reversed(half_span), strict=True)]
    return f'<path d="{polyline(inner + outer)} Z"/>'


def rungs(motion: str) -> str:
    out = []
    half = RUNG_PITCH * RUNG_DUTY / 2
    for m in range(CROSSINGS):
        centre = (m + 0.5) * DT  # middle of a lobe, where the strands are widest
        for j in range(RUNGS_PER_LOBE):
            t = centre + (j - (RUNGS_PER_LOBE - 1) / 2) * RUNG_PITCH
            if 2 * AMP * abs(math.sin(K * t)) - STROKE >= RUNG_MIN_VISIBLE:
                out.append(rung(t - half, t + half))
    return f'<g fill="{RUNG}">{"".join(out)}{motion}</g>'


def fan(u0: float, u1: float) -> str:
    """Wedges coloured by polar angle over [u0, u1], covering the thread's band."""
    r_in, r_out = R - AMP - STROKE, R + AMP + STROKE
    n = math.ceil((u1 - u0) / WEDGE)
    out = []
    for i in range(n):
        a = u0 + i * WEDGE
        b = a + WEDGE * (1 + WEDGE_OVERLAP)
        pts = [point_at(a, r_in), point_at(b, r_in), point_at(b, r_out), point_at(a, r_out)]
        out.append(f'<path d="{polyline(pts)} Z" fill="{colour(a + WEDGE / 2)}"/>')
    return "".join(out)


def segment(index: int, sweep: float, motion: str) -> tuple[str, str]:
    """Mask and paint for the thread between two consecutive under-passes.

    Passes alternate over/under along the thread: even passes (t = m * DT)
    go over, odd ones under. Splitting at under-passes puts every joint
    inside a cut, so segments meet without a visible seam. The fan covers
    every angle the mask reaches while it turns by `sweep`."""
    m0 = 2 * index + 1
    t0, t1 = m0 * DT, (m0 + 2) * DT
    reach = CUT_REACH * DT
    # The other pass through the same crossing is CROSSINGS passes away
    # (half the thread), and CROSSINGS is odd, so it is an over-pass.
    cuts = "".join(
        f'<path d="{thread(tc - reach, tc + reach, 48)}"/>'
        for tc in (((m + CROSSINGS) % PASSES) * DT for m in (m0, m0 + 2))
    )
    mask = (
        f'<mask id="ow-seg{index}" maskUnits="userSpaceOnUse" '
        f'x="0" y="0" width="{fmt(SIZE)}" height="{fmt(SIZE)}">'
        f'<g fill="none">'
        f'<path d="{thread(t0, t1, 160)}" stroke="#FFF" stroke-width="{fmt(STROKE)}"/>'
        f'<g stroke="#000" stroke-width="{fmt(STROKE + 2 * GAP)}">{cuts}</g>'
        f"{motion}</g></mask>"
    )
    paint = f'<g mask="url(#ow-seg{index})">{fan(t0 - FAN_PAD, t1 + sweep + FAN_PAD)}</g>'
    return mask, paint


def build(animated: bool) -> str:
    sweep = SPIN if animated else 0.0
    spin = rotation(SPIN_PERIOD, math.degrees(SPIN)) if animated else ""
    turn = rotation(TURN_PERIOD, 360.0) if animated else ""
    masks, paints = zip(*(segment(i, sweep, spin) for i in range(CROSSINGS)), strict=True)
    size = fmt(SIZE)
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}" '
            f'width="{size}" height="{size}" role="img" aria-label="OmicsWeave">',
            "<title>OmicsWeave</title>",
            f"<defs>{''.join(masks)}</defs>",
            "<g>",
            rungs(spin),
            *paints,
            f"{turn}</g>",
            "</svg>\n",
        ]
    )


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Draw the OmicsWeave logo.")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=root / "assets",
        help="directory for logo.svg and logo-animated.svg (default: assets/)",
    )
    args = parser.parse_args()
    if CROSSINGS % 2 == 0 or not math.isclose(2 * K, CROSSINGS):
        parser.error("K must be a half-integer (one thread, odd crossing count)")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for name, animated in (("logo.svg", False), ("logo-animated.svg", True)):
        path = args.out_dir / name
        path.write_text(build(animated))
        print(path)


if __name__ == "__main__":
    main()
