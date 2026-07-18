#!/usr/bin/python3
"""Generate Supplementary Figure S2: transmembrane-gateway method.

The figure is intentionally methodological. Empirical gateway distributions
remain in main Figure 4, while this supplement defines how gateway records are
constructed and aggregated.
"""
from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
from pathlib import Path

from reportlab.graphics import renderPDF, renderSVG
from reportlab.graphics.shapes import Circle, Drawing, Line, Polygon, Rect, String
from reportlab.lib.colors import HexColor, white
from reportlab.pdfbase import pdfmetrics


W, H = 510, 360  # 7.08 x 5.00 inches
INK = HexColor("#20242B")
MUTED = HexColor("#69727D")
GRID = HexColor("#DDE3E8")
PALE = HexColor("#F5F7F9")
BLUE = HexColor("#4C72B0")
BLUE_PALE = HexColor("#EAF0F8")
TEAL = HexColor("#2A9D8F")
TEAL_PALE = HexColor("#E8F5F2")
ORANGE = HexColor("#E07A3F")
ORANGE_PALE = HexColor("#FCF0E8")
RED = HexColor("#C44E52")
GREEN = HexColor("#55A868")
GREEN_PALE = HexColor("#EAF4EC")
VIOLET = HexColor("#8172B2")
VIOLET_PALE = HexColor("#F0EDF7")
GREY_PALE = HexColor("#EEF1F3")

FONT_REGULAR = "Helvetica"
FONT_BOLD = "Helvetica-Bold"
TEXT_SCALE = 1.08


class Canvas:
    """Small top-left-coordinate wrapper around ReportLab vector shapes."""

    def __init__(self) -> None:
        self.d = Drawing(W, H)
        self.d.add(Rect(0, 0, W, H, fillColor=white, strokeColor=None))

    @staticmethod
    def fy(y: float) -> float:
        return H - y

    def line(self, x1, y1, x2, y2, color=INK, width=0.8, dash=None) -> None:
        shape = Line(x1, self.fy(y1), x2, self.fy(y2), strokeColor=color, strokeWidth=width)
        if dash:
            shape.strokeDashArray = dash
        self.d.add(shape)

    def rect(self, x, y, w, h, fill=white, stroke=GRID, width=0.7, radius=4) -> None:
        self.d.add(
            Rect(
                x,
                self.fy(y + h),
                w,
                h,
                rx=radius,
                ry=radius,
                fillColor=fill,
                strokeColor=stroke,
                strokeWidth=width,
            )
        )

    def circle(self, x, y, r, fill=white, stroke=INK, width=0.8) -> None:
        self.d.add(
            Circle(x, self.fy(y), r, fillColor=fill, strokeColor=stroke, strokeWidth=width)
        )

    def text(self, x, y, value, size=7, color=INK, bold=False, anchor="start") -> None:
        self.d.add(
            String(
                x,
                self.fy(y),
                value,
                fontName=FONT_BOLD if bold else FONT_REGULAR,
                fontSize=size * TEXT_SCALE,
                fillColor=color,
                textAnchor=anchor,
            )
        )

    def lines(self, x, y, values, size=6.4, color=INK, bold=False, leading=8, anchor="start") -> None:
        for i, value in enumerate(values):
            self.text(x, y + i * leading, value, size, color, bold, anchor)

    def arrow(self, x1, y1, x2, y2, color=INK, width=0.8, head=3.0) -> None:
        self.line(x1, y1, x2, y2, color, width)
        angle = math.atan2(y2 - y1, x2 - x1)
        p1 = (x2, y2)
        p2 = (
            x2 - head * math.cos(angle) + head * 0.55 * math.sin(angle),
            y2 - head * math.sin(angle) - head * 0.55 * math.cos(angle),
        )
        p3 = (
            x2 - head * math.cos(angle) - head * 0.55 * math.sin(angle),
            y2 - head * math.sin(angle) + head * 0.55 * math.cos(angle),
        )
        points = []
        for x, y in (p1, p2, p3):
            points.extend((x, self.fy(y)))
        self.d.add(Polygon(points, fillColor=color, strokeColor=color, strokeWidth=0.2))

    def double_arrow(self, x1, y1, x2, y2, color=INK, width=0.7, head=2.4) -> None:
        self.arrow(x1, y1, x2, y2, color, width, head)
        self.arrow(x2, y2, x1, y1, color, width, head)


def panel_heading(c: Canvas, x: float, label: str, title: str, subtitle: str) -> None:
    c.text(x, 20, label, 11.5, INK, True)
    c.text(x + 17, 20, title, 8.8, INK, True)
    c.text(x + 17, 33, subtitle, 6.2, MUTED)


def pmax_text(
    c: Canvas,
    x: float,
    y: float,
    size: float,
    color=INK,
    bold: bool = False,
    anchor: str = "start",
    prefix: str = "",
    suffix: str = "",
) -> None:
    """Draw p with an upright 'max' subscript and optional surrounding text."""
    font = FONT_BOLD if bold else FONT_REGULAR
    sub_size = size * 0.68
    draw_size = size * TEXT_SCALE
    draw_sub_size = sub_size * TEXT_SCALE
    prefix_w = pdfmetrics.stringWidth(prefix, font, draw_size)
    p_w = pdfmetrics.stringWidth("p", font, draw_size)
    sub_w = pdfmetrics.stringWidth("max", font, draw_sub_size)
    suffix_w = pdfmetrics.stringWidth(suffix, font, draw_size)
    total_w = prefix_w + p_w + sub_w + suffix_w + 0.5
    start = x
    if anchor == "middle":
        start -= total_w / 2
    elif anchor == "end":
        start -= total_w
    c.text(start, y, prefix, size, color, bold)
    p_x = start + prefix_w
    c.text(p_x, y, "p", size, color, bold)
    c.text(p_x + p_w, y + draw_size * 0.23, "max", sub_size, color, bold)
    c.text(p_x + p_w + sub_w + 0.5, y, suffix, size, color, bold)


def draw_panel_a(c: Canvas) -> None:
    x0, x1 = 7, 166
    panel_heading(c, x0, "A", "Seven TM-helix gateways", "Top view; TM assignments from GPCRdb")

    cx, cy, radius, hr = 84, 127, 48, 10.5
    points = []
    for i in range(7):
        angle = math.radians(-90 + i * 360 / 7)
        points.append((cx + radius * math.cos(angle), cy + radius * math.sin(angle)))

    # Gateway edges sit behind the helices and establish the seven-record ring.
    for i in range(7):
        p1, p2 = points[i], points[(i + 1) % 7]
        c.line(p1[0], p1[1], p2[0], p2[1], TEAL, 4.0)

    c.circle(cx, cy, 22, BLUE_PALE, GRID, 0.7)
    c.lines(cx, cy - 3, ["receptor", "core"], 6.2, INK, True, 7.2, "middle")

    for i, (hx, hy) in enumerate(points, start=1):
        c.circle(hx, hy, hr, white, BLUE, 2.1)
        c.text(hx, hy + 2.2, f"TM{i}", 6.3, BLUE, True, "middle")

    gateway_labels = [
        (130, 69, "TM1-TM2"),
        (150, 113, "TM2-TM3"),
        (136, 172, "TM3-TM4"),
        (84, 194, "TM4-TM5"),
        (31, 172, "TM5-TM6"),
        (17, 113, "TM6-TM7"),
        (38, 69, "TM7-TM1"),
    ]
    for lx, ly, label in gateway_labels:
        c.text(lx, ly, label, 5.4, TEAL, True, "middle")

    c.rect(x0 + 5, 214, x1 - x0 - 10, 40, BLUE_PALE, None, 0, 4)
    c.text(x0 + 13, 228, "One record per neighboring helix pair", 6.4, INK, True)
    c.text(x0 + 13, 240, "7 interfaces x 4 reported metrics", 6.5, MUTED)

    c.text(x0 + 12, 275, "Identity contract", 6.5, BLUE, True)
    c.line(x0 + 12, 281, x1 - 12, 281, GRID, 0.7)
    c.text(x0 + 12, 294, "Same interface keys in every included system", 6.2, INK)
    c.text(x0 + 12, 307, "Missing record: null + reason code", 6.2, INK)
    c.text(x0 + 12, 320, "Never substitute a numerical zero", 6.2, RED, True)


def draw_panel_b(c: Canvas) -> None:
    x0, x1 = 174, 355
    panel_heading(c, x0, "B", "Frame-level classification", "Two decisions are applied to each sampled frame")

    c.text(x0 + 5, 55, "1", 8.2, white, True, "middle")
    c.circle(x0 + 5, 52.5, 7, TEAL, TEAL, 0.5)
    c.text(x0 + 17, 55, "WEDGED LIPID ATOM", 7.0, TEAL, True)

    # Side view: dual-helix proximity plus the central z-band.
    c.rect(x0 + 13, 75, 150, 46, GREEN_PALE, None, 0, 3)
    c.text(x0 + 88, 82, "central TM z-band", 5.7, GREEN, True, "middle")
    c.text(x0 + 88, 117, "25th-75th percentile", 5.4, GREEN, anchor="middle")
    c.rect(x0 + 24, 61, 22, 81, white, BLUE, 1.5, 10)
    c.rect(x0 + 132, 61, 22, 81, white, BLUE, 1.5, 10)
    c.text(x0 + 35, 151, "helix i", 5.6, BLUE, True, "middle")
    c.text(x0 + 143, 151, "helix i+1", 5.6, BLUE, True, "middle")
    atom_x, atom_y = x0 + 90, 100
    c.circle(atom_x, atom_y, 4.2, ORANGE, white, 1.0)
    c.double_arrow(x0 + 46, atom_y, atom_x - 5, atom_y, TEAL, 0.8, 2.2)
    c.double_arrow(atom_x + 5, atom_y, x0 + 132, atom_y, TEAL, 0.8, 2.2)
    c.text(x0 + 66, atom_y - 4, "<5.0 \u00c5", 5.2, TEAL, True, "middle")
    c.text(x0 + 111, atom_y - 4, "<5.0 \u00c5", 5.2, TEAL, True, "middle")
    c.circle(x0 + 90, 60, 3.5, GREY_PALE, MUTED, 0.7)
    c.text(x0 + 96, 61.5, "outside z-band", 5.2, MUTED)

    c.rect(x0 + 17, 160, 147, 24, TEAL_PALE, None, 0, 3)
    c.text(x0 + 90, 170, "distance to BOTH helices < 5.0 \u00c5", 6.0, INK, True, "middle")
    c.text(x0 + 90, 179, "AND atom z lies inside the central band", 5.9, INK, True, "middle")

    c.line(x0 + 8, 195, x1 - 7, 195, GRID, 0.8)
    c.text(x0 + 5, 214, "2", 8.2, white, True, "middle")
    c.circle(x0 + 5, 211.5, 7, ORANGE, ORANGE, 0.5)
    c.text(x0 + 17, 214, "OPEN FRAME", 7.0, ORANGE, True)

    # Top view: penetration is the deepest inward displacement in this frame.
    wall_x, center_x = x0 + 147, x0 + 85
    c.line(wall_x, 229, wall_x, 286, ORANGE, 1.3, [3, 2])
    c.text(wall_x - 2, 224, "nearest helix-pair C-alpha wall", 5.3, ORANGE, True, "end")
    c.circle(center_x, 258, 8, BLUE_PALE, BLUE, 1.0)
    c.text(center_x, 260, "core", 4.8, BLUE, True, "middle")
    for ax, ay in ((x0 + 132, 245), (x0 + 119, 268), (x0 + 102, 255)):
        c.circle(ax, ay, 3.4, GREY_PALE, MUTED, 0.6)
    deepest_x, deepest_y = x0 + 102, 255
    c.circle(deepest_x, deepest_y, 4.1, RED, white, 0.9)
    c.double_arrow(deepest_x + 5, deepest_y, wall_x - 1, deepest_y, RED, 1.0, 2.5)
    pmax_text(c, (deepest_x + wall_x) / 2, deepest_y - 5, 5.6, RED, True, "middle")
    c.text(x0 + 17, 239, "wedged atoms", 5.4, MUTED)
    c.arrow(x0 + 61, 237, x0 + 116, 246, MUTED, 0.5, 2.0)
    pmax_text(
        c,
        x0 + 90,
        287,
        5.4,
        MUTED,
        anchor="middle",
        suffix=" = deepest inward displacement",
    )

    c.rect(x0 + 17, 296, 147, 28, ORANGE_PALE, None, 0, 3)
    pmax_text(c, x0 + 90, 307, 7.0, RED, True, "middle", suffix=" >= 0.5 \u00c5")
    c.text(x0 + 90, 317, "frame classified as OPEN", 6.2, INK, True, "middle")


def draw_replica_lane(c: Canvas, y: float, label: str, color, open_pattern) -> None:
    x0 = 365
    c.text(x0 + 7, y + 3, label, 6.2, color, True)
    start = x0 + 42
    for i, is_open in enumerate(open_pattern):
        fill = RED if is_open else GREY_PALE
        stroke = RED if is_open else MUTED
        c.rect(start + i * 7.0, y - 5, 5, 10, fill, stroke, 0.45, 1)
    c.arrow(start + 46, y, x0 + 99, y, MUTED, 0.6, 2.2)
    c.rect(x0 + 101, y - 8, 34, 16, color if color != VIOLET else VIOLET, None, 0, 3)
    c.text(x0 + 118, y + 2, "summary", 5.2, white, True, "middle")


def metric_box(c: Canvas, x: float, y: float, title: str, detail: str, fill) -> None:
    c.rect(x, y, 64, 30, fill, None, 0, 3)
    c.text(x + 5, y + 11, title, 5.7, INK, True)
    c.text(x + 5, y + 22, detail, 5.0, MUTED)


def draw_panel_c(c: Canvas) -> None:
    x0, x1 = 363, 505
    panel_heading(c, x0, "C", "Replica aggregation", "Frames -> 3 replicas -> system record")

    c.text(x0 + 7, 55, "sampled frames", 5.6, MUTED)
    c.text(x0 + 119, 55, "4 metrics", 5.6, MUTED, anchor="middle")
    draw_replica_lane(c, 75, "rep 1", BLUE, [0, 1, 0, 0, 1, 1])
    draw_replica_lane(c, 106, "rep 2", GREEN, [0, 0, 1, 0, 1, 0])
    draw_replica_lane(c, 137, "rep 3", VIOLET, [1, 0, 1, 1, 0, 0])
    c.text(x0 + 52, 153, "closed", 5.1, MUTED)
    c.rect(x0 + 42, 148, 5, 7, GREY_PALE, MUTED, 0.4, 1)
    c.text(x0 + 93, 153, "open", 5.1, MUTED)
    c.rect(x0 + 83, 148, 5, 7, RED, RED, 0.4, 1)

    bus_x = x1 - 2
    for start_y in (75, 106, 137):
        c.line(x0 + 135, start_y, bus_x, start_y, GRID, 0.7)
    c.line(bus_x, 75, bus_x, 172, GRID, 0.7)
    c.line(bus_x, 172, x0 + 70, 172, GRID, 0.7)
    c.arrow(x0 + 70, 172, x0 + 70, 185, MUTED, 0.8, 2.5)

    c.rect(x0 + 10, 187, 124, 36, BLUE_PALE, BLUE, 0.8, 4)
    c.text(x0 + 72, 199, "SYSTEM-LEVEL GATEWAY RECORD", 6.2, BLUE, True, "middle")
    c.text(x0 + 72, 211, "mean of the 3 replica summaries", 6.0, INK, True, "middle")
    c.text(x0 + 72, 220, "for each interface and metric", 5.4, MUTED, anchor="middle")

    c.arrow(x0 + 72, 224, x0 + 72, 235, MUTED, 0.8, 2.5)
    c.rect(x0 + 19, 237, 106, 27, VIOLET_PALE, VIOLET, 0.8, 4)
    c.text(x0 + 72, 248, "95% bootstrap interval", 6.2, VIOLET, True, "middle")
    c.text(x0 + 72, 258, "1,000 resamples; fixed seed", 5.4, MUTED, anchor="middle")

    c.rect(x0 + 3, 274, 64, 30, ORANGE_PALE, None, 0, 3)
    c.text(x0 + 8, 285, "penetration", 5.7, INK, True)
    pmax_text(c, x0 + 8, 296, 5.0, MUTED, prefix="mean ", suffix=" (\u00c5)")
    metric_box(c, x0 + 72, 274, "penetration_p90", "90th percentile (\u00c5)", ORANGE_PALE)
    metric_box(c, x0 + 3, 309, "open_fraction", "open frames / all frames", TEAL_PALE)
    metric_box(c, x0 + 72, 309, "occupancy", "wedged atoms / frame", GREEN_PALE)


def build_figure() -> Drawing:
    c = Canvas()
    c.line(169.5, 10, 169.5, 340, GRID, 0.8)
    c.line(358.5, 10, 358.5, 340, GRID, 0.8)
    draw_panel_a(c)
    draw_panel_b(c)
    draw_panel_c(c)
    return c.d


def validate_source_data(data_dir: Path) -> None:
    path = data_dir / "Supplementary_Figure_S2_gateway_method_source.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 12
    values = {row["parameter_or_metric"]: row["value_or_definition"] for row in rows}
    expected = {
        "helix-pair distance cutoff": "5.0",
        "central TM z-band lower bound": "25",
        "central TM z-band upper bound": "75",
        "open-frame penetration threshold": "0.5",
        "adjacent TM interfaces": "7",
        "source replicas": "3",
        "bootstrap confidence level": "95",
        "bootstrap resamples": "1000",
        "penetration": "mean p_max",
        "penetration_p90": "90th percentile",
        "open_fraction": "open frames / all frames",
        "occupancy": "wedged atoms / frame",
    }
    assert values == expected


def main() -> None:
    repository_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=repository_root / "source_data",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository_root / "figures",
    )
    args = parser.parse_args()
    validate_source_data(args.data_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure = build_figure()
    stem = args.output_dir / "supplementary_figure_s2_gateway_method"
    renderPDF.drawToFile(figure, str(stem.with_suffix(".pdf")))
    renderSVG.drawToFile(figure, str(stem.with_suffix(".svg")))
    subprocess.run(
        [
            shutil.which("pdftoppm") or "/usr/bin/pdftoppm",
            "-r",
            "600",
            "-png",
            "-singlefile",
            str(stem.with_suffix(".pdf")),
            str(stem),
        ],
        check=True,
    )
    for path in sorted(args.output_dir.glob("supplementary_figure_s2_gateway_method.*")):
        print(path)


if __name__ == "__main__":
    main()
