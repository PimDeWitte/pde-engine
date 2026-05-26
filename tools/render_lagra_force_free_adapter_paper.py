#!/usr/bin/env python3
"""Render the Lagra/pde-engine force-free adapter paper.

The visual structure mirrors LagraPaperV1: pipeline, vocabulary chips,
graph anatomy, per-step data flow, exploded kernel, verifier triangle, and
derivation-style measurement check.  The content is intentionally narrow:
it documents the force-free adapter and its measurements, not a full
pde-engine reproduction.
"""

from __future__ import annotations

import argparse
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
HTML_OUT = DOCS / "lagra-force-free-adapter-paper.html"
PDF_OUT = DOCS / "lagra-force-free-adapter-paper.pdf"


GREEN = "#139f73"
TEAL = "#15839c"
PURPLE = "#6c48c9"
BLUE = "#3f4fc6"
MAGENTA = "#c6246b"
ORANGE = "#c97606"
INK = "#243043"
MUTED = "#697386"
PAPER = "#fbfaf7"
LINE = "#9aa5b1"


def rect(
    x: int,
    y: int,
    w: int,
    h: int,
    fill: str,
    stroke: str | None = None,
    rx: int = 6,
    cls: str = "",
) -> str:
    stroke_attr = f' stroke="{stroke}"' if stroke else ""
    cls_attr = f' class="{cls}"' if cls else ""
    return (
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" '
        f'fill="{fill}"{stroke_attr}{cls_attr}/>'
    )


def text(
    x: int,
    y: int,
    s: str,
    size: int = 13,
    fill: str = INK,
    anchor: str = "middle",
    weight: str = "400",
    style: str = "",
    family: str = "Georgia, serif",
) -> str:
    style_attr = f' style="{style}"' if style else ""
    return (
        f'<text x="{x}" y="{y}" font-size="{size}" fill="{fill}" '
        f'text-anchor="{anchor}" font-weight="{weight}" '
        f'font-family="{family}"{style_attr}>{s}</text>'
    )


def arrow(x1: int, y1: int, x2: int, y2: int, color: str = INK, dash: bool = False) -> str:
    dash_attr = ' stroke-dasharray="5 5"' if dash else ""
    return (
        f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" '
        f'stroke="{color}" stroke-width="1.8"{dash_attr} marker-end="url(#arrow)"/>'
    )


def defs() -> str:
    return f"""
    <defs>
      <marker id="arrow" viewBox="0 0 10 10" refX="8.5" refY="5"
        markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 0 L 10 5 L 0 10 z" fill="{INK}"/>
      </marker>
      <marker id="arrowMagenta" viewBox="0 0 10 10" refX="8.5" refY="5"
        markerWidth="6" markerHeight="6" orient="auto-start-reverse">
        <path d="M 0 0 L 10 5 L 0 10 z" fill="{MAGENTA}"/>
      </marker>
    </defs>
    """


def chip(x: int, y: int, w: int, h: int, color: str, title: str, lines: list[str]) -> str:
    out = [rect(x, y, w, h, color, rx=5)]
    out.append(text(x + w // 2, y + 22, title, 13, "white", weight="700"))
    for i, line in enumerate(lines):
        out.append(text(x + w // 2, y + 42 + 18 * i, line, 11, "white", family="Courier New, monospace"))
    return "\n".join(out)


def figure(num: int, title: str, svg: str, caption: str) -> str:
    figure_class = ' class="pagebreak"' if num == 5 else ""
    return f"""
    <figure{figure_class}>
      {svg}
      <figcaption><strong>Figure {num}: {title}.</strong> {caption}</figcaption>
    </figure>
    """


def svg1_pipeline() -> str:
    boxes = [
        (34, 88, 106, 58, GREEN, ".lagra", "adapter scene"),
        (174, 88, 118, 58, TEAL, "Lagrangian", "Graph"),
        (326, 88, 126, 58, PURPLE, "Action", "Adapter"),
        (486, 88, 126, 58, BLUE, "pde-engine", "Validator"),
        (646, 88, 110, 58, MAGENTA, "Verifier", "gates"),
    ]
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 235" role="img" aria-label="Lagra adapter pipeline">',
        defs(),
        text(210, 54, "compile-time", 12, MUTED),
        text(552, 54, "per-run measurement loop", 12, MUTED),
    ]
    for x, y, w, h, color, a, b in boxes:
        parts.append(rect(x, y, w, h, color, rx=5))
        parts.append(text(x + w // 2, y + 25, a, 14, "white", weight="700"))
        parts.append(text(x + w // 2, y + 44, b, 12, "white", weight="700"))
    parts.extend(
        [
            arrow(140, 117, 174, 117),
            arrow(292, 117, 326, 117),
            arrow(452, 117, 486, 117),
            arrow(612, 117, 646, 117),
            text(87, 168, "typed proof-search", 12, INK),
            text(87, 184, "constants and gates", 12, INK),
            text(233, 168, "vertices = checks", 12, INK),
            text(233, 184, "edges = measurements", 12, INK),
            text(389, 168, "imports symbols,", 12, INK),
            text(389, 184, "builds det M", 12, INK),
            text(549, 168, "valid / invalid,", 12, INK),
            text(549, 184, "cache hashes", 12, INK),
            text(701, 168, "PASS -> paper", 12, INK),
            text(701, 184, "FAIL -> boundary", 12, INK),
            '<path d="M 704 146 L 704 207 L 389 207 L 389 147" fill="none" '
            f'stroke="{MAGENTA}" stroke-width="1.4" stroke-dasharray="5 5" '
            'marker-end="url(#arrowMagenta)"/>',
            text(546, 226, "paper output is admitted only after integrity + readiness gates pass", 12, MAGENTA),
            "</svg>",
        ]
    )
    return "\n".join(parts)


def svg2_vocab() -> str:
    return "\n".join(
        [
            '<svg class="figure-svg" viewBox="0 0 790 270" role="img" aria-label="Adapter vocabulary chips">',
            defs(),
            text(394, 24, "AdapterKind = RegistrySolution | NegativeControl | PointCheck", 12, INK, family="Courier New, monospace"),
            text(394, 43, "| RotationContext | HealthSmoke | PaperArtifact", 12, INK, family="Courier New, monospace"),
            chip(54, 62, 205, 72, GREEN, "RegistrySolution", ["input: KNOWN_SOLUTIONS", "measure: det M == 0"]),
            chip(292, 62, 205, 72, PURPLE, "NegativeControl", ["input: rho*z, exp(rho*z)", "measure: det M != 0"]),
            chip(531, 62, 205, 72, ORANGE, "PointCheck", ["input: rational points", "measure: abs(det M)"]),
            chip(54, 158, 205, 72, TEAL, "RotationContext", ["input: Omega in cache key", "measure: valid/invalid split"]),
            chip(292, 158, 205, 72, BLUE, "HealthSmoke", ["input: bounded CLI run", "measure: rc0 / counts"]),
            chip(531, 158, 205, 72, MAGENTA, "PaperArtifact", ["input: JSON catalog", "measure: exact hashes"]),
            text(394, 252, "Every chip has a typed input, a deterministic measurement, and an admitted output field.", 12, MUTED),
            "</svg>",
        ]
    )


def svg3_graph() -> str:
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 355" role="img" aria-label="Anatomy of Lagra adapter graph">',
        defs(),
        text(154, 52, "external inputs", 12, MUTED, style="font-style:italic"),
        text(625, 52, "measured outputs", 12, MUTED, style="font-style:italic"),
        rect(48, 78, 152, 46, GREEN, rx=5),
        text(124, 105, "pde-engine registry", 12, "white", weight="700"),
        rect(48, 147, 152, 46, PURPLE, rx=5),
        text(124, 174, "control expressions", 12, "white", weight="700"),
        rect(48, 216, 152, 46, TEAL, rx=5),
        text(124, 243, "context: Omega", 12, "white", weight="700"),
        rect(310, 82, 170, 46, PAPER, stroke=LINE, rx=5),
        text(395, 110, "symbols rho, z", 12, INK),
        rect(310, 151, 170, 46, PAPER, stroke=LINE, rx=5),
        text(395, 179, "A, B, L_T operators", 12, INK),
        rect(310, 220, 170, 46, PAPER, stroke=LINE, rx=5),
        text(395, 248, "determinant M", 12, INK),
        arrow(200, 101, 310, 105),
        arrow(200, 170, 310, 174),
        arrow(200, 239, 310, 243),
        arrow(395, 128, 395, 151),
        arrow(395, 197, 395, 220),
        rect(594, 73, 154, 50, GREEN, rx=5),
        text(671, 102, "7/7 exact zero", 12, "white", weight="700"),
        rect(594, 139, 154, 50, PURPLE, rx=5),
        text(671, 168, "2/2 exact nonzero", 12, "white", weight="700"),
        rect(594, 205, 154, 50, TEAL, rx=5),
        text(671, 234, "Omega cache split", 12, "white", weight="700"),
        rect(594, 271, 154, 50, MAGENTA, rx=5),
        text(671, 300, "paper/integrity pass", 12, "white", weight="700"),
        arrow(480, 105, 594, 98),
        arrow(480, 174, 594, 164),
        arrow(480, 243, 594, 230),
        arrow(480, 243, 594, 296),
        '<path d="M 268 308 H 522" fill="none" stroke="#333" stroke-width="1.2"/>',
        text(395, 300, "Lagra edge sum:", 12, MUTED),
        text(395, 326, "E = E_registry + E_controls + E_rotation + E_health + E_artifact", 13, INK, family="Courier New, monospace"),
        "</svg>",
    ]
    return "\n".join(parts)


def svg4_dataflow() -> str:
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 270" role="img" aria-label="One adapter run as data flow">',
        defs(),
        rect(40, 62, 130, 48, TEAL, rx=5),
        text(105, 83, "Lagra catalog", 12, "white", weight="700"),
        text(105, 101, "constants", 11, "white"),
        rect(40, 145, 130, 48, ORANGE, rx=5),
        text(105, 166, "Current input", 12, "white", weight="700"),
        text(105, 184, "u(rho,z), Omega", 11, "white"),
        rect(244, 100, 130, 60, PURPLE, rx=5),
        text(309, 122, "Exact adapter", 12, "white", weight="700"),
        text(309, 141, "u -> det M", 11, "white"),
        rect(436, 100, 132, 60, PAPER, stroke=LINE, rx=5),
        text(502, 122, "Measurements", 12, INK, weight="700"),
        text(502, 141, "zero/nonzero/abs", 11, INK),
        rect(632, 100, 126, 60, BLUE, rx=5),
        text(695, 122, "Output JSON", 12, "white", weight="700"),
        text(695, 141, "same type each run", 11, "white"),
        arrow(170, 86, 244, 119),
        arrow(170, 169, 244, 141),
        arrow(374, 130, 436, 130),
        arrow(568, 130, 632, 130),
        text(107, 128, "explicit inputs", 11, MUTED),
        text(309, 187, "rebuilds det M for this context", 11, MUTED),
        text(502, 187, "exact + point measurements", 11, MUTED),
        text(695, 187, "readiness-audit input", 11, MUTED),
        "</svg>",
    ]
    return "\n".join(parts)


def svg5_kernel() -> str:
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 355" role="img" aria-label="Adapter kernel exploded">',
        defs(),
        rect(42, 62, 142, 50, TEAL, rx=5),
        text(113, 84, "Gate graph", 12, "white", weight="700"),
        text(113, 102, "active checks", 11, "white"),
        rect(42, 139, 142, 50, ORANGE, rx=5),
        text(113, 161, "State", 12, "white", weight="700"),
        text(113, 179, "pde checkout + data", 11, "white"),
        rect(260, 69, 310, 116, "#f4f0ff", stroke=PURPLE, rx=6),
        text(282, 94, "Adapter kernel", 13, PURPLE, "start", weight="700"),
        text(282, 120, "for check in graph.active_checks():", 12, INK, "start", family="Courier New, monospace"),
        text(302, 140, "measurement = check.run(context)", 12, INK, "start", family="Courier New, monospace"),
        text(302, 160, "artifact.accumulate(measurement)", 12, INK, "start", family="Courier New, monospace"),
        rect(634, 101, 110, 52, MAGENTA, rx=5),
        text(689, 123, "Artifact", 12, "white", weight="700"),
        text(689, 141, "metrics", 11, "white"),
        arrow(184, 87, 260, 105),
        arrow(184, 164, 260, 149),
        arrow(570, 127, 634, 127),
        text(604, 113, "accumulate", 11, MAGENTA, style="font-style:italic"),
        text(395, 218, "each check.run(context) dispatches to one of:", 12, MUTED),
        chip(42, 236, 160, 70, GREEN, "Known", ["det M == 0", "count -> 7"]),
        chip(221, 236, 160, 70, PURPLE, "Control", ["det M != 0", "count -> 2"]),
        chip(400, 236, 160, 70, TEAL, "Cache", ["hash0 != hash1", "valid != invalid"]),
        chip(579, 236, 160, 70, BLUE, "Health", ["112 generated", "rc0, no Lean fail"]),
        text(394, 333, "Adding a new adapter measurement is adding a check, not editing the loop.", 12, INK),
        "</svg>",
    ]
    return "\n".join(parts)


def svg6_triangle() -> str:
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 330" role="img" aria-label="Verifier triangle for Lagra and pde-engine">',
        defs(),
        rect(70, 40, 310, 250, "#f7fbf8", rx=6),
        rect(410, 40, 310, 250, "#fbf9f5", rx=6),
        text(225, 66, "(a) Lagra adapter verifier", 13, INK, weight="700"),
        text(565, 66, "(b) Full pde-engine reproduction", 13, INK, weight="700"),
        rect(95, 88, 260, 48, "#e8fbf2", stroke=GREEN, rx=5),
        text(114, 110, "proposition", 12, GREEN, "start", weight="700"),
        text(114, 128, "7/7 det M zero; cache split holds", 11, INK, "start", family="Courier New, monospace"),
        rect(115, 162, 220, 48, "#f4efff", stroke=PURPLE, rx=5),
        text(134, 184, "witness", 12, PURPLE, "start", weight="700"),
        text(134, 202, "JSON metrics + source hashes", 11, INK, "start"),
        rect(115, 235, 220, 38, "#fff0f6", stroke=MAGENTA, rx=5),
        text(134, 259, "checker: integrity gate -> PASS", 11, INK, "start", weight="700"),
        arrow(225, 136, 225, 162),
        arrow(225, 210, 225, 235),
        rect(435, 88, 260, 48, "#e8fbf2", stroke=GREEN, rx=5),
        text(454, 110, "proposition", 12, GREEN, "start", weight="700"),
        text(454, 128, "all validators + Lean path reproduce", 11, INK, "start", family="Courier New, monospace"),
        rect(455, 162, 220, 48, "#f4efff", stroke=PURPLE, rx=5),
        text(474, 184, "witness", 12, PURPLE, "start", weight="700"),
        text(474, 202, "bounded smoke: rc0, 2 known", 11, INK, "start"),
        rect(455, 235, 220, 38, "#fff0f6", stroke=MAGENTA, rx=5),
        text(474, 259, "checker: full 7-solution path not green", 11, INK, "start", weight="700"),
        arrow(565, 136, 565, 162),
        arrow(565, 210, 565, 235),
        text(394, 313, "Same discipline as the Lagra verifier triangle; different claim admitted.", 12, MUTED),
        "</svg>",
    ]
    return "\n".join(parts)


def svg7_derivation() -> str:
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 360" role="img" aria-label="Derivation-style cache regression">',
        defs(),
        rect(250, 34, 290, 56, "#f4efff", stroke=PURPLE, rx=5),
        text(395, 58, "u = rho^2 exp(-2 z)", 14, INK, weight="700"),
        text(395, 78, "same expression, two validation contexts", 12, INK),
        text(152, 115, "Omega = 0 branch", 12, MUTED, style="font-style:italic"),
        rect(62, 132, 190, 52, "#e8fbf2", stroke=GREEN, rx=5),
        text(157, 154, "pde-engine validator", 12, INK, weight="700"),
        text(157, 172, "valid == True", 11, INK),
        rect(62, 216, 190, 52, "#e8fbf2", stroke=GREEN, rx=5),
        text(157, 238, "cache hash H0", 12, INK, weight="700"),
        text(157, 256, "context includes Omega=0", 11, INK),
        text(638, 115, "Omega = 1 branch", 12, MUTED, style="font-style:italic"),
        rect(538, 132, 190, 52, "#fff0f6", stroke=MAGENTA, rx=5),
        text(633, 154, "pde-engine validator", 12, INK, weight="700"),
        text(633, 172, "valid == False", 11, INK),
        rect(538, 216, 190, 52, "#fff0f6", stroke=MAGENTA, rx=5),
        text(633, 238, "cache hash H1", 12, INK, weight="700"),
        text(633, 256, "context includes Omega=1", 11, INK),
        arrow(325, 90, 175, 132),
        arrow(465, 90, 615, 132),
        arrow(157, 184, 157, 216),
        arrow(633, 184, 633, 216),
        rect(256, 294, 278, 46, "#dff7ef", stroke=GREEN, rx=5),
        text(395, 316, "measured invariant: H0 != H1 and rotating det M != 0", 12, INK, weight="700"),
        text(395, 334, "det M = 2048 rho^9 exp(-12 z)", 12, INK, family="Courier New, monospace"),
        arrow(252, 242, 316, 294),
        arrow(538, 242, 474, 294),
        "</svg>",
    ]
    return "\n".join(parts)


def svg8_discoveries() -> str:
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 360" role="img" aria-label="Novel discovery artifact flow">',
        defs(),
        rect(44, 50, 162, 58, TEAL, rx=5),
        text(125, 73, "Bounded engine", 12, "white", weight="700"),
        text(125, 91, "depth 2, 112 rows", 11, "white"),
        rect(314, 50, 162, 58, PURPLE, rx=5),
        text(395, 73, "Filter", 12, "white", weight="700"),
        text(395, 91, "valid, non-known", 11, "white"),
        rect(584, 50, 162, 58, GREEN, rx=5),
        text(665, 73, "Independent check", 12, "white", weight="700"),
        text(665, 91, "rebuild det M", 11, "white"),
        arrow(206, 79, 314, 79),
        arrow(476, 79, 584, 79),
        text(260, 65, "76 valid", 11, MUTED),
        text(530, 65, "both rho,z", 11, MUTED),
        chip(54, 150, 205, 70, BLUE, "rho + z", ["oblique linear", "det M -> 0"]),
        chip(292, 150, 205, 70, MAGENTA, "rho/z", ["self-similar", "det M -> 0"]),
        chip(531, 150, 205, 70, ORANGE, "rho^2 + z^2", ["quadratic radius", "det M -> 0"]),
        chip(54, 248, 205, 70, GREEN, "rho^2 + z", ["parabolic polynomial", "det M -> 0"]),
        chip(292, 248, 205, 70, TEAL, "-rho^2 + z^2 + 1", ["hyperbolic quadratic", "det M -> 0"]),
        chip(531, 248, 205, 70, PURPLE, "rho/(1 - z)", ["rational sheet", "det M -> 0"]),
        text(394, 340, "Novelty is scoped to the repo registry: not identical to the seven registered known_solutions.", 12, INK),
        "</svg>",
    ]
    return "\n".join(parts)


def svg9_kerr_gate() -> str:
    parts = [
        '<svg class="figure-svg" viewBox="0 0 790 380" role="img" aria-label="Kerr paper target gate">',
        defs(),
        rect(38, 48, 158, 60, BLUE, rx=5),
        text(117, 72, "Literature target", 12, "white", weight="700"),
        text(117, 91, "Kerr FFE / GSE", 11, "white"),
        rect(232, 48, 158, 60, PURPLE, rx=5),
        text(311, 72, "Strict gate", 12, "white", weight="700"),
        text(311, 91, "paper criteria", 11, "white"),
        rect(426, 48, 158, 60, TEAL, rx=5),
        text(505, 72, "Current engine", 12, "white", weight="700"),
        text(505, 91, "linear surrogate", 11, "white"),
        rect(620, 48, 132, 60, MAGENTA, rx=5),
        text(686, 72, "Artifact", 12, "white", weight="700"),
        text(686, 91, "no_candidate_yet", 11, "white"),
        arrow(196, 78, 232, 78),
        arrow(390, 78, 426, 78),
        arrow(584, 78, 620, 78),
        rect(56, 150, 170, 64, "#eef3ff", stroke=BLUE, rx=5),
        text(76, 174, "paper target", 12, BLUE, "start", weight="700"),
        text(76, 194, "stationary axisymmetric", 11, INK, "start"),
        text(76, 210, "magnetically dominated", 11, INK, "start"),
        rect(310, 136, 176, 112, "#f4efff", stroke=PURPLE, rx=5),
        text(330, 160, "admission checks", 12, PURPLE, "start", weight="700"),
        text(330, 182, "depends on r, x, a", 11, INK, "start"),
        text(330, 198, "finite denominators", 11, INK, "start"),
        text(330, 214, "small-spin anchor", 11, INK, "start"),
        text(330, 230, "axis / horizon regular", 11, INK, "start"),
        rect(552, 150, 176, 100, "#fff0f6", stroke=MAGENTA, rx=5),
        text(572, 174, "measured output", 12, MAGENTA, "start", weight="700"),
        text(572, 194, "306 rows + 568 corrections", 11, INK, "start"),
        text(572, 210, "48-coeff solve: EmptySet", 11, INK, "start"),
        text(572, 226, "BZ O(a^2) anchor: pass", 11, INK, "start"),
        text(572, 242, "878 inputs, 0 admitted", 11, INK, "start"),
        arrow(226, 182, 310, 192),
        arrow(486, 192, 552, 188),
        rect(248, 306, 294, 46, "#fbf9f5", stroke=LINE, rx=5),
        text(395, 328, "guardrail: surrogate rows cannot become Kerr paper claims", 12, INK, weight="700"),
        text(395, 346, "full residual exists; no expression passes it", 11, MUTED),
        arrow(505, 108, 430, 306, dash=True),
        "</svg>",
    ]
    return "\n".join(parts)


def html() -> str:
    figs = [
        figure(
            1,
            "The Lagra/pde-engine adapter pipeline in one figure",
            svg1_pipeline(),
            "This is the direct analogue of the Lagra pipeline figure. The input is a typed Lagra proof-search source plus a local pde-engine checkout. The measured output is not a trajectory; it is a JSON evidence artifact and an admitted paper section. The dashed return path is the strict gate: failed measurements do not become prose.",
        ),
        figure(
            2,
            "The adapter vocabulary",
            svg2_vocab(),
            "The original Lagra vocabulary figure lists LagrangianNode variants. Here the closed vocabulary is the set of adapter checks. Each variant names exactly what enters the check and exactly what is measured on the way out.",
        ),
        figure(
            3,
            "Anatomy of the Lagra adapter graph",
            svg3_graph(),
            "The external pde-engine registry, control expressions, and validation context enter as graph-visible inputs. The edges are measurement obligations. The outputs are the counts and facts consumed by the proof-search catalog: exact-zero solutions, exact-nonzero controls, rotation-cache separation, and paper-readiness status.",
        ),
        figure(
            4,
            "One adapter run as a data-flow diagram",
            svg4_dataflow(),
            "This mirrors the one-timestep figure in the Lagra paper. The action is not a mechanical time step; it is one deterministic validation pass. The adapter rebuilds the symbolic determinant for the current expression and context, measures exact and pointwise quantities, and writes a typed JSON artifact.",
        ),
        figure(
            5,
            "The adapter kernel, exploded",
            svg5_kernel(),
            "The loop is deliberately boring. Heterogeneity lives in check.run(context), not in the loop. Adding a Kerr no-promotion guard, a force-free cache regression, or a future PDE boundary check means adding a new typed check with a measurement contract.",
        ),
        figure(
            6,
            "Verifier triangle for the adapter and the full reproduction boundary",
            svg6_triangle(),
            "This copies the structural shape of the Lagra verifier triangle while keeping the trust event honest. The adapter proposition passes with JSON witnesses and source hashes. The bounded pde-engine smoke is now process-clean: return code 0, no Lean build failure, 112 generated expressions, 76 valid rows, and 2 known vertical forms. The full seven-solution reproduction proposition still remains non-green because the bounded depth-2 smoke does not recover all seven known solutions.",
        ),
        figure(
            7,
            "Verifiability by derivation for the cache regression",
            svg7_derivation(),
            "The same expression is sent down two independent validation-context branches. Lagra measures the branch split, the pde-engine validator result, and the exact rotating determinant. The conclusion is a code obligation: the validator cache key must include Omega and validation mode.",
        ),
        figure(
            8,
            "Engine-discovered candidate foliations admitted by the adapter",
            svg8_discoveries(),
            "The bounded depth-2 pde-engine run produces valid rows that are not identical to the seven registered known solutions. The exporter does not trust the cache: it rebuilds the determinant independently and admits only both-variable candidates whose determinant simplifies exactly to zero. This is an engine-novel result relative to the repository registry, not a literature-priority claim.",
        ),
        figure(
            9,
            "Kerr paper-target gate",
            svg9_kerr_gate(),
            "The next target uses the same measurement discipline in the negative direction. The bounded Kerr run generated 306 rows; the strict gate scans those rows, 456 targeted finite-spin corrections, 112 Kerr-metric-resummed corrections, and four probes against the closed split-monopole nonlinear Kerr Grad-Shafranov residual. It also solves the 48-coefficient leading-order ansatz and gets EmptySet, verifies the known O(a^2) Blandford-Znajek perturbative anchor, then emits no_candidate_yet because no finite-spin exact expression passes the criteria.",
        ),
    ]

    figures = "\n".join(figs)
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Lagra force-free adapter for pde-engine</title>
  <style>
    @page {{ size: A4; margin: 17mm 18mm; }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: #f0eee9;
      color: {INK};
      font-family: Georgia, "Times New Roman", serif;
      font-size: 13.5px;
      line-height: 1.44;
    }}
    main {{
      width: 210mm;
      margin: 0 auto;
      background: #fffefa;
      padding: 17mm 18mm 20mm;
      min-height: 297mm;
    }}
    .running {{ color: #a0a0a0; font-size: 10px; margin-bottom: 28px; }}
    h1 {{
      text-align: center;
      font-size: 18px;
      line-height: 1.25;
      margin: 0 0 24px;
    }}
    .byline {{ text-align: center; font-size: 12px; margin-bottom: 36px; }}
    .abstract-title {{ text-align: center; font-weight: 700; margin: 0 0 12px; }}
    .abstract {{ max-width: 620px; margin: 0 auto 30px; text-align: justify; }}
    h2 {{ font-size: 16px; margin: 26px 0 8px; }}
    h3 {{ font-size: 14px; margin: 20px 0 6px; }}
    p {{ text-align: justify; margin: 0 0 10px; }}
    code, pre {{ font-family: "Courier New", monospace; }}
    pre {{
      background: #f3f2ed;
      padding: 10px 12px;
      overflow-wrap: anywhere;
      white-space: pre-wrap;
      font-size: 11.5px;
      line-height: 1.35;
    }}
    figure {{ margin: 25px 0 26px; break-inside: avoid; }}
    .pagebreak {{ break-before: page; page-break-before: always; }}
    .figure-svg {{
      width: 100%;
      height: auto;
      display: block;
      margin: 0 auto 8px;
    }}
    figcaption {{ font-size: 12.5px; line-height: 1.35; text-align: left; }}
    .math {{ font-family: Georgia, "Times New Roman", serif; font-style: italic; }}
    table {{ width: 100%; border-collapse: collapse; margin: 14px 0; font-size: 12px; }}
    th, td {{ border-bottom: 1px solid #ddd8ce; padding: 6px 7px; text-align: left; vertical-align: top; }}
    th {{ background: #f6f2e8; }}
    .claim {{ border-left: 3px solid {MAGENTA}; padding-left: 10px; }}
    @media print {{
      body {{ background: #fffefa; }}
      main {{ width: auto; margin: 0; padding: 0; }}
      a {{ color: inherit; text-decoration: none; }}
    }}
  </style>
</head>
<body>
<main>
  <div class="running">Lagra force-free adapter for pde-engine</div>
  <h1>Lagra as a measurement harness for pde-engine force-free validation</h1>
  <div class="byline">Pim de Witte · Lagra Project · pde-engine adapter note · May 26, 2026</div>

  <div class="abstract-title">Abstract</div>
  <p class="abstract">
    This note rewrites the pde-engine force-free adapter in the visual style of
    <em>Lagra: A differentiable physics kernel for verifiable domain operations</em>.
    The purpose is narrower than the Lagra kernel paper.  We do not claim a new
    literature-priority result or a full reproduction of the pde-engine discovery
    pipeline.  We show how Lagra enters a local pde-engine checkout, which
    inputs are admitted into the Lagra proof-search graph, which quantities are
    measured, and which claims are allowed into the paper.  The positive result
    is exact and useful: seven pde-engine registry force-free foliations simplify
    to determinant zero; two negative controls simplify to nonzero determinants;
    three rational checkpoints separate known solutions from controls; and a
    rotation-context cache regression proves that the validator cache key must
    include <code>Omega</code> and the validation mode; and a bounded engine pass
    exports six non-registered, independently rechecked candidate foliations
    as a process-control baseline rather than as a publication target.
    A second gate records the next Kerr paper target, appends targeted one-term
    and two-term finite-spin correction grammars, runs a 48-coefficient
    leading-order solve, checks the literature slow-rotation perturbative
    anchor, and emits <code>no_candidate_yet</code> rather than promoting the
    current linear surrogate.
    The negative result is equally
    important: the bounded pde-engine smoke is now process-clean, but the full
    seven-solution pde-engine/Lean reproduction path remains non-green on this
    checkout and is recorded as a boundary, not upgraded into a proof.
  </p>

  {figures}

  <h2>1. What exactly enters Lagra</h2>
  <p>
    The adapter has four external inputs.  First, it imports the pde-engine
    force-free symbols <code>rho</code> and <code>z</code> from
    <code>problems/force_free/validator.py</code>.  Second, it loads the seven known
    Compere force-free foliations from the pde-engine problem registry.  Third,
    it supplies two false controls, <code>rho*z</code> and <code>exp(rho*z)</code>.  Fourth,
    it runs one rotation-context check on <code>rho**2*exp(-2*z)</code> with
    <code>Omega=0</code> and <code>Omega=1</code> against a shared cache database.
  </p>
  <p>
    Nothing enters through prose.  Each input is bound into a Lagra proof-search
    artifact and then checked by <code>publication_readiness_audit.py</code> and
    <code>proof_search_integrity_gate.py</code>.  The paper is downstream of those
    gates.  This is the central point of using Lagra here: the artifact catalog,
    not the narrative, decides whether the claim is admissible.
  </p>

  <h2>2. What exactly is measured</h2>
  <p>
    For an expression <code>u(rho,z)</code>, the adapter evaluates the force-free
    determinant from Compere, Gralla, and Lupsasca, <em>Force-Free
    Foliations</em>, Phys. Rev. D 94, 124012 (2016), Eq. 2.14 / Section 2.4:
  </p>
  <pre>det([[L_T(A), L_T(B)], [L_T^2(A), L_T^2(B)]])</pre>
  <pre>A = u_rho_rho + u_z_z - u_rho/rho
B = u_rho**2 + u_z**2
T = u_z*d_rho - u_rho*d_z</pre>
  <p>
    The strict measurements are:
  </p>
  <table>
    <tr><th>Measurement</th><th>Observed result</th><th>Claim admitted</th></tr>
    <tr><td>Registry exact determinant</td><td><code>7/7</code> simplify to <code>0</code></td><td>The Lagra adapter matches the pde-engine registry solutions for this determinant.</td></tr>
    <tr><td>Negative-control exact determinant</td><td><code>rho*z -> 16*rho*z</code>; <code>exp(rho*z) -> 16*rho*z*exp(6*rho*z)</code></td><td>The operator is not vacuously zero.</td></tr>
    <tr><td>Rational point checks</td><td>Known solutions below <code>1e-60</code>; controls above <code>1e-6</code></td><td>The exact result survives independent numeric checkpoints.</td></tr>
    <tr><td>Rotation-context cache regression</td><td><code>rho**2*exp(-2*z)</code> valid at <code>Omega=0</code>, invalid at <code>Omega=1</code></td><td>The pde-engine cache key must include rotation context.</td></tr>
    <tr><td>Rotating determinant</td><td><code>2048*rho**9*exp(-12*z)</code></td><td>The rotating failure is a symbolic fact, not only a cache artifact.</td></tr>
  </table>

  <h2>3. Engine process-control candidates admitted by the measurement</h2>
  <p>
    The PR now includes a second artifact, generated by
    <code>tools/export_force_free_novel_discoveries.py</code>.  It runs the bounded
    pde-engine search, mines rows that are valid but not registered as one of
    the seven known force-free solutions, and then independently rebuilds the
    determinant instead of trusting the validator cache.  The claim is scoped
    narrowly: these are process-control candidates, engine-novel relative to the
    repository registry, not claimed as literature-first discoveries.
  </p>
  <table>
    <tr><th>Engine-discovered expression</th><th>Class</th><th>Independent witness</th></tr>
    <tr><td><code>rho + z</code></td><td>linear two-coordinate foliation</td><td><code>det M = 0</code></td></tr>
    <tr><td><code>rho/z</code></td><td>scale-invariant ratio foliation</td><td><code>det M = 0</code></td></tr>
    <tr><td><code>rho**2 + z**2</code></td><td>quadratic radius foliation</td><td><code>det M = 0</code></td></tr>
    <tr><td><code>rho**2 + z</code></td><td>parabolic polynomial foliation</td><td><code>det M = 0</code></td></tr>
    <tr><td><code>-rho**2 + z**2 + 1</code></td><td>hyperbolic quadratic foliation</td><td><code>det M = 0</code></td></tr>
    <tr><td><code>rho/(1 - z)</code></td><td>rational geometric-sum foliation</td><td><code>det M = 0</code></td></tr>
  </table>
  <p>
    The machine-readable version is
    <code>docs/force-free-novel-discoveries.json</code>; the review-friendly
    version is <code>docs/force-free-novel-discoveries.md</code>.
  </p>
  <p>
    The full process record and next research target are documented in
    <code>docs/discovery-process-and-next-target.md</code>.  The real paper target is
    not this force-free registry problem; it is a Kerr force-free
    Grad-Shafranov/extreme-Kerr target where the literature does not already
    provide an exact analytic solution meeting the same criteria.
  </p>

  <h2>4. Kerr paper-target gate</h2>
  <p>
    The next target gate is now a concrete artifact, not only a paragraph:
    <code>docs/kerr-paper-target-gate.json</code> and
    <code>docs/kerr-paper-target-gate.md</code>.  It ran the current
    <code>kerr_magnetosphere</code> engine harness at depth 2, recorded
    <code>306</code> generated rows and <code>306</code> completed validations,
    then scanned all <code>306</code> generated expressions, <code>192</code>
    one-term correction candidates, <code>264</code> two-term correction
    candidates, and four probes against the
    full gate.  It admitted <code>0</code> paper candidates.
  </p>
  <p>
    This is the right failure mode.  The current repository target is only a
    linear surrogate for generation, while paper admission is now decided by
    the closed split-monopole nonlinear residual from Mahlmann et al. Eq.
    <code>GSLightCylinder</code>, rewritten with <code>x = cos(theta)</code>.  The
    strict gate requires dependence on <code>r</code>, <code>x</code>, and
    <code>a</code>; finite denominators on rational safe points; a small-spin
    anchor to <code>1 - x</code> or <code>x</code>; rejection of those exact anchors
    as known; and exact full-residual zero after the point checks.  No row meets
    all of those criteria, so the artifact status is <code>no_candidate_yet</code>.
  </p>
  <p>
    The targeted correction grammar is deliberately bounded:
    <code>Psi = 1 - x + a**2*c*basis(r,x)</code>, with
    <code>c in {-1, -1/2, 1/2, 1}</code>, eight angular factors, and six radial
    factors.  That Cartesian product generates <code>192</code> correction
    candidates.  The two-term screen uses a 12-function sub-basis and sign
    coefficients to generate <code>264</code> more candidates.  All
    <code>456</code> finite-spin correction candidates pass the cheap strict
    prechecks before the full residual, and all <code>456</code> fail exact-zero
    residual.  The follow-up metric-resummed screen adds <code>112</code>
    Kerr-denominator candidates using <code>Sigma</code>, <code>Delta</code>,
    and the light-cylinder metric block; all <code>112</code> pass strict
    prechecks and all <code>112</code> fail exact-zero residual.
  </p>
  <p>
    The stricter screen is not a larger fixed-coefficient search.  It solves the
    leading <code>a**2</code> residual coefficient for the full one-term basis:
    <code>Psi = 1 - x + a**2*sum_i c_i*basis_i(r,x)</code> with <code>M = 1</code>.
    That produces <code>66</code> polynomial equations in <code>48</code> unknown
    coefficients; the matrix has shape <code>[66, 48]</code>, and
    <code>linsolve</code> returns <code>EmptySet</code>.  So even before the full
    nonlinear residual, this ansatz family has no leading-order correction.
  </p>
  <p>
    The next calibration is positive but deliberately perturbative.  The gate
    encodes the Blandford-Znajek split-monopole slow-rotation correction from
    Tanabe-Nagataki and Pan-Yu:
    <code>Psi = 1 - x + a**2*x*(1-x**2)*R(r)</code>, where <code>R(r)</code>
    contains the literature logarithm and dilogarithm radial terms.  After
    applying <code>polylog(1,z) = -log(1-z)</code>, the coefficient of
    <code>a**2</code> in the closed residual simplifies exactly to <code>0</code>.
    This proves the gate recognizes the known <code>O(a**2)</code> perturbative
    anchor, while still keeping it out of the finite-spin exact candidate set.
  </p>
  <p>
    The JSON artifact now includes a criteria-status matrix.  It marks the
    literature target, bounded engine generation, strict prechecks, nonlinear
    residual screen, finite-spin correction screens, leading-order coefficient
    solve, and literature perturbative anchor calibration as implemented.  It
    deliberately marks broader equivalence filters and global
    horizon/axis/light-surface regularity as not sufficient for a positive
    paper claim.  That distinction matters: this is a target packet and
    negative gate for a no-known-exact-solution problem, not a positive exact
    Kerr solution.
  </p>
  <p>
    The literature motivation is explicitly separated from the current
    surrogate.  Mahlmann et al. frame static, axisymmetric, force-free Kerr
    magnetospheres around the relativistic Grad-Shafranov equation and numerical
    solution methods.  Camilloni et al. state that for extreme Kerr there is no
    known exact analytic stationary, axisymmetric, magnetically dominated
    force-free solution.  The gate is designed for that claim family, not for
    more examples from the Compere force-free foliation registry.
  </p>

  <h2>5. The pde-engine patch forced by the measurement</h2>
  <p>
    The old validator cache hashed only the expression string.  That was wrong:
    a non-rotating result could be replayed in a rotating context.  The patch
    therefore includes the expression, schema, <code>Omega</code>,
    <code>check_regularity</code>, <code>fast_point_only</code>, and <code>use_lean</code> in the hash.
    The patch also removes the old fast-point placeholder derivative shortcut.
    Fast point mode still avoids the expensive full path, but it now builds the
    exact symbolic determinant before stopping at the exact point check.
  </p>
  <p>
    The parallel reproduction worker had a separate bug: spawned children tried
    to import <code>physics_agent.problems</code>, which is not the package namespace in
    this repository.  The patch switches the child import to
    <code>problems.load_problem</code> and adds a local fallback import for
    <code>PreciseFoliationValidator</code>.
  </p>

  <h2>6. Reproduction commands</h2>
  <pre>python3 -m py_compile \
  general_method_paper_reproduction.py \
  problems/force_free/validator.py \
  tests/test_force_free_validator_cache.py \
  tools/export_force_free_novel_discoveries.py \
  tools/export_kerr_paper_target_gate.py \
  tests/test_force_free_novel_discoveries.py \
  tests/test_kerr_paper_target_gate.py

python3 -m unittest \
  tests/test_force_free_validator_cache.py \
  tests/test_force_free_novel_discoveries.py \
  tests/test_kerr_paper_target_gate.py

python3 tools/export_force_free_novel_discoveries.py \
  --run-engine --max-depth 2 --validators 1 \
  --timeout-s 90 --validation-timeout-s 3

python3 tools/export_kerr_paper_target_gate.py \
  --run-engine --max-depth 2 --validators 1 \
  --timeout-s 90 --validation-timeout-s 3</pre>
  <p>
    The Lagra-side commands that produced the measurements were:
  </p>
  <pre>python3 experiments/proof-search/pde_engine_force_free_point_gate.py \
  --pde-engine-root /Users/p/dev/pde-engine

python3 experiments/proof-search/pde_engine_reproduction_health_gate.py \
  --pde-engine-root /Users/p/dev/pde-engine \
  --timeout-seconds 20

python3 experiments/proof-search/render_yukawa_finite_range_paper.py
python3 experiments/proof-search/publication_readiness_audit.py
python3 experiments/proof-search/proof_search_integrity_gate.py</pre>

  <h2>7. References</h2>
  <p>
    [1] Geoffrey Compere, Samuel E. Gralla, Alexandru Lupsasca,
    <em>Force-Free Foliations</em>, Phys. Rev. D 94, 124012 (2016),
    arXiv:1606.06727, DOI:10.1103/PhysRevD.94.124012.
    The paper formulates force-free electrodynamics in terms of field-line
    foliations; in the stationary axisymmetric case the object used here is a
    foliation of the half-plane.
  </p>
  <p>
    Links: <code>https://arxiv.org/abs/1606.06727</code> and
    <code>https://doi.org/10.1103/PhysRevD.94.124012</code>.
  </p>
  <p>
    [2] J. F. Mahlmann, P. Cerda-Duran, M. A. Aloy et al.,
    <em>Numerically solving the relativistic Grad-Shafranov equation in Kerr
    spacetimes: Numerical techniques</em>, MNRAS 477, 3927-3946 (2018),
    arXiv:1802.00815, DOI:10.1093/mnras/sty858.
  </p>
  <p>
    [3] F. Camilloni, G. Grignani, T. Harmark, R. Oliveri, M. Orselli,
    <em>Moving away from the Near-Horizon Attractor of the Extreme Kerr
    Force-Free Magnetosphere</em>, JCAP 10, 048 (2020), arXiv:2007.15665,
    DOI:10.1088/1475-7516/2020/10/048.
  </p>
  <p>
    [4] K. Tanabe and S. Nagataki,
    <em>Higher Order Terms of Kerr Parameter for Blandford-Znajek Monopole
    Solution</em>, arXiv:0802.0908.
  </p>
  <p>
    [5] Z. Pan and C. Yu,
    <em>Fourth-order split monopole perturbation solutions to the
    Blandford-Znajek mechanism</em>, arXiv:1503.05248.
  </p>

  <h2>8. Boundary</h2>
  <p class="claim">
    The adapter is a positive exact-symbolic and pointwise Lagra measurement of
    the pde-engine force-free boundary.  It is not a Lean theorem.  It is not a
    full pde-engine paper reproduction.  The bounded depth-2 pde-engine smoke is
    now process-clean, and it yields six independently rechecked process-control
    candidates not identical to the registered known-solution functions.  The
    full reproduction boundary remains: the bounded pass only finds two known
    vertical canonical forms and does not recover the seven known Compere
    solutions.  The Kerr gate is a separate <code>no_candidate_yet</code> artifact:
    it prevents the current linear surrogate from being promoted into the next
    paper target, appends targeted finite-spin corrections, and checks the
    closed split-monopole nonlinear residual directly.  Its leading-order
    coefficient solve also returns <code>EmptySet</code> for the 48-function
    one-term ansatz.  Its literature anchor screen does recognize the known
    <code>O(a**2)</code> Blandford-Znajek perturbative correction.
  </p>
</main>
</body>
</html>
"""


def render_pdf(html_path: Path, pdf_path: Path) -> None:
    try:
        from playwright.sync_api import sync_playwright
    except Exception as exc:  # pragma: no cover - environment-dependent
        raise RuntimeError(
            "Playwright is required for PDF rendering; HTML was still written."
        ) from exc

    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_page(viewport={"width": 1240, "height": 1754})
        page.goto(html_path.resolve().as_uri(), wait_until="networkidle")
        page.pdf(path=str(pdf_path), format="A4", print_background=True)
        browser.close()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-pdf", action="store_true", help="Only write the HTML artifact.")
    args = parser.parse_args()

    DOCS.mkdir(parents=True, exist_ok=True)
    content = "\n".join(line.rstrip() for line in html().splitlines()) + "\n"
    HTML_OUT.write_text(content, encoding="utf-8")
    if not args.no_pdf:
        render_pdf(HTML_OUT, PDF_OUT)
    print(HTML_OUT)
    if PDF_OUT.exists():
        print(PDF_OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
