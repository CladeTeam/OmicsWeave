"""Draw the OmicsWeave logo and write it to assets/logo.svg.

The mark is a single thread that circles twice before it closes, so its two
laps read as the complementary strands of a double helix bent into a ring
(the torus knot T(2, 2K) for half-integer K). Its colour cycles through the
four base colours once over both laps, so the two laps at any angle sit half
a cycle apart: A faces T and C faces G.

Over/under is drawn by removing the lower pass wherever it comes within GAP
of the upper pass. Nothing is painted in a background colour, so the rungs
stay visible and the logo has a transparent background.

Usage: python scripts/plot_logo.py [-o OUT]
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

# The colour gradient is painted as short overlapping pieces; each segment's
# mask then cuts one clean outline, so piece edges never show.
PIECE_STEPS = 72  # pieces per segment
PIECE_OVERLAP = 0.35  # fraction of a step each piece extends into the next
PIECE_WIDTH = STROKE + 6.0
CUT_REACH = 0.45  # cut length either side of a crossing, as a fraction of DT

CROSSINGS = round(2 * K)
PASSES = 2 * CROSSINGS  # each crossing is visited twice along the thread
DT = math.pi / K  # parameter step between consecutive crossing passes
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


def rung(t0: float, t1: float, steps: int = 4) -> str:
    """Sector between angles t0 and t1, spanning the two strands.

    Its ends follow the strand centre lines, so they stay hidden under the
    strands drawn on top."""
    ts = [t0 + (t1 - t0) * i / steps for i in range(steps + 1)]
    half_span = [AMP * abs(math.sin(K * t)) - RUNG_TUCK for t in ts]
    inner = [point_at(t, R - s) for t, s in zip(ts, half_span, strict=True)]
    outer = [point_at(t, R + s) for t, s in zip(reversed(ts), reversed(half_span), strict=True)]
    return f'<path d="{polyline(inner + outer)} Z"/>'


def rungs() -> str:
    out = []
    half = RUNG_PITCH * RUNG_DUTY / 2
    for m in range(CROSSINGS):
        centre = (m + 0.5) * DT  # middle of a lobe, where the strands are widest
        for j in range(RUNGS_PER_LOBE):
            t = centre + (j - (RUNGS_PER_LOBE - 1) / 2) * RUNG_PITCH
            if 2 * AMP * abs(math.sin(K * t)) - STROKE >= RUNG_MIN_VISIBLE:
                out.append(rung(t - half, t + half))
    return f'<g fill="{RUNG}">{"".join(out)}</g>'


def segment(index: int) -> tuple[str, str]:
    """Mask and paint for the thread between two consecutive under-passes.

    Passes alternate over/under along the thread: even passes (t = m * DT)
    go over, odd ones under. Splitting at under-passes puts every joint
    inside a cut, so segments meet without a visible seam."""
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
        f'x="0" y="0" width="{fmt(SIZE)}" height="{fmt(SIZE)}" fill="none">'
        f'<path d="{thread(t0, t1, 160)}" stroke="#FFF" stroke-width="{fmt(STROKE)}"/>'
        f'<g stroke="#000" stroke-width="{fmt(STROKE + 2 * GAP)}">{cuts}</g>'
        f"</mask>"
    )

    step = (t1 - t0) / PIECE_STEPS
    pieces = "".join(
        f'<path d="{thread(a, a + step * (1 + PIECE_OVERLAP), 2)}" '
        f'stroke="{colour(a + step / 2)}"/>'
        for a in (t0 + i * step for i in range(PIECE_STEPS))
    )
    paint = (
        f'<g mask="url(#ow-seg{index})" fill="none" stroke-width="{fmt(PIECE_WIDTH)}">{pieces}</g>'
    )
    return mask, paint


def build() -> str:
    masks, paints = zip(*(segment(i) for i in range(CROSSINGS)), strict=True)
    size = fmt(SIZE)
    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {size} {size}" '
            f'width="{size}" height="{size}" role="img" aria-label="OmicsWeave">',
            "<title>OmicsWeave</title>",
            f"<defs>{''.join(masks)}</defs>",
            rungs(),
            *paints,
            "</svg>\n",
        ]
    )


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description="Draw the OmicsWeave logo.")
    parser.add_argument(
        "-o",
        "--out",
        type=Path,
        default=root / "assets" / "logo.svg",
        help="output SVG path (default: assets/logo.svg)",
    )
    args = parser.parse_args()
    if CROSSINGS % 2 == 0 or not math.isclose(2 * K, CROSSINGS):
        parser.error("K must be a half-integer (one thread, odd crossing count)")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(build())
    print(args.out)


if __name__ == "__main__":
    main()
