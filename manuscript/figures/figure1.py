"""
Figure 1 mockup: Complementary frameworks for synthetic scenario design.

Tool choice: matplotlib (v3.10.7, confirmed available).
  - Exports SVG + PDF natively; both needed for manuscript workflow.
  - Code is the regenerable authoritative source — spec changes become
    one-line edits rather than manual SVG surgery.
  - drawsvg/SVG authoring would be faster for a pure one-off hand-drawn
    schematic, but cannot be regenerated reproducibly when the spec changes,
    which is likely during review iterations.
  - matplotlib's FancyBboxPatch + annotate + figure-coordinate placement
    gives adequate quality for a methods schematic at column width.

Outputs:  manuscript/figures/figure1_mockup.svg
          manuscript/figures/figure1_mockup.pdf
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
from pathlib import Path

# ── Palette ──────────────────────────────────────────────────────────────────
BLUE   = '#4472C4'   # input / forcing space
ORANGE = '#ED7D31'   # generator / simulation model
GREEN  = '#548235'   # hazard / outcome space
PURPLE = '#8064A2'   # analysis / scenario discovery
GRAY   = '#B8B8B8'   # ensemble (light: black text readable on it)
RED    = '#C00000'   # MOEA optimizer (panel D only)
GOLD   = '#C8940A'   # DESIGN label (darkened for contrast on white bg)
LGRAY  = '#F3F3F3'   # panel background
SLGN   = '#EAF1E6'   # strip background (pale green tint)
WHITE  = '#FFFFFF'
DKGRAY = '#1A1A1A'

# ── Figure setup ─────────────────────────────────────────────────────────────
FW, FH = 14.0, 7.4
fig = plt.figure(figsize=(FW, FH), facecolor='white', dpi=150)
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, FW)
ax.set_ylim(0, FH)
ax.axis('off')

# ── Layout constants (all in figure-inches) ───────────────────────────────────
LM       = 0.20          # left margin
RM       = 0.20          # right margin
PGAP     = 0.13          # gap between panels
PW       = (FW - LM - RM - 3*PGAP) / 4   # ≈ 3.3225 in
BW       = PW - 0.34     # standard box width (0.17 in margin each side)
BH       = 0.57          # standard box height
R_PAD    = 0.03          # rounded-corner pad for FancyBboxPatch

STRIP_BOT = 0.10
STRIP_H   = 0.78
PAN_BOT   = STRIP_BOT + STRIP_H + 0.16   # ≈ 1.04
PAN_TOP   = FH - 0.08                     # ≈ 7.32
PH        = PAN_TOP - PAN_BOT             # ≈ 6.28

# Panel left-x edges
PLX = [LM + i*(PW + PGAP) for i in range(4)]
PCX = [pl + PW/2 for pl in PLX]          # panel center x

# Box y-centers from top down (shared across panels A-C; panel D deviates)
TITLE_Y = PAN_TOP - 0.50
Y = [
    PAN_TOP - 1.15,   # Y[0] first box
    PAN_TOP - 2.05,   # Y[1]
    PAN_TOP - 2.95,   # Y[2]
    PAN_TOP - 3.85,   # Y[3]
    PAN_TOP - 4.75,   # Y[4] last box
]
# Gap between bottom of box[i] and top of box[i+1] ≈ 0.9 - 0.57 = 0.33 in — room for arrow label

# ── Low-level drawing primitives ─────────────────────────────────────────────

def rbox(cx, cy, w, h, fc, ec='white', lw=0.8, zorder=3, alpha=1.0):
    """Draw centered rounded rectangle."""
    p = FancyBboxPatch(
        (cx - w/2, cy - h/2), w, h,
        boxstyle=f'round,pad={R_PAD}',
        facecolor=fc, edgecolor=ec, linewidth=lw,
        zorder=zorder, alpha=alpha,
    )
    ax.add_patch(p)

def txt(x, y, s, fs=6.8, fc=WHITE, fw='bold', ha='center', va='center', z=4,
        style='normal', wrap=False):
    ax.text(x, y, s, ha=ha, va=va, fontsize=fs, color=fc,
            fontweight=fw, zorder=z, style=style,
            multialignment='center')

def box_and_text(cx, cy, w, h, fc, main_txt, sub_txt=None,
                 main_fs=6.6, sub_fs=5.5, text_color=WHITE, ec='white', lw=0.8):
    """Draw a box with a main label and optional italic sub-label."""
    rbox(cx, cy, w, h, fc, ec=ec, lw=lw)
    if sub_txt:
        txt(cx, cy + 0.10, main_txt, fs=main_fs, fc=text_color)
        txt(cx, cy - 0.13, sub_txt, fs=sub_fs, fc=text_color, fw='normal', style='italic')
    else:
        txt(cx, cy, main_txt, fs=main_fs, fc=text_color)

def arrow(x1, y1, x2, y2, color='black', lw=1.2, dashed=False,
          rad=0.0, lbl=None, lbl_side='right', lbl_fs=5.7, lbl_color=None, z=5):
    """Annotate an arrow; optionally add a labelled badge alongside."""
    ls = (0, (4, 2.5)) if dashed else 'solid'
    ax.annotate(
        '', xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle='->', color=color, lw=lw,
            linestyle=ls,
            connectionstyle=f'arc3,rad={rad}',
            shrinkA=2, shrinkB=2,
        ),
        zorder=z,
    )
    if lbl:
        # Badge position: midpoint, offset to avoid overlap with arrow shaft
        mx = (x1 + x2) / 2
        my = (y1 + y2) / 2
        x_off = +0.52 if lbl_side == 'right' else -0.52
        lc = lbl_color or GOLD
        ax.text(mx + x_off, my, lbl,
                ha='center', va='center', fontsize=lbl_fs,
                color=lc, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.12', fc=WHITE, ec=lc,
                          lw=0.8, alpha=0.97),
                zorder=8)

def panel_bg_and_title(pi, title, letter):
    # Background
    ax.add_patch(plt.Rectangle(
        (PLX[pi], PAN_BOT), PW, PH,
        fc=LGRAY, ec='#C0C0C0', lw=0.8, zorder=1))
    # Panel letter badge (top-left)
    ax.text(PLX[pi] + 0.11, PAN_TOP - 0.14, letter,
            ha='center', va='center', fontsize=8.5,
            fontweight='bold', color=DKGRAY, zorder=5)
    # Title
    ax.text(PLX[pi] + PW/2, TITLE_Y, title,
            ha='center', va='center', fontsize=7.5,
            fontweight='bold', color=DKGRAY,
            multialignment='center', zorder=5)

# ── Panel backgrounds & titles ────────────────────────────────────────────────
TITLES = [
    '(a)  Forward exploratory\nmodeling',
    '(b)  Bottom-up vulnerability\nassessment',
    '(c)  Library-and-subsample',
    '(d)  MOEA-FIND:\nhazard-space targeting',
]
LETTERS = ['a', 'b', 'c', 'd']
for pi in range(4):
    panel_bg_and_title(pi, TITLES[pi], LETTERS[pi])

# ── PANEL A — Forward exploratory modeling ────────────────────────────────────
PI = 0
cx = PCX[PI]

box_and_text(cx, Y[0], BW, BH, BLUE,
             'Uncertain parameters /\nHistorical forcing',
             sub_txt='uncertainty space  X')
box_and_text(cx, Y[1], BW, BH, ORANGE,
             'Streamflow generator /\nclimate model', None)
box_and_text(cx, Y[2], BW, BH+0.04, GRAY,
             'Large SOW ensemble',
             sub_txt='(many realizations)', text_color=DKGRAY, ec='#909090', lw=0.6)
box_and_text(cx, Y[3], BW, BH, ORANGE,
             'Simulation model', sub_txt='(Pywr-DRB)')
box_and_text(cx, Y[4], BW, BH, GREEN,
             'Performance metrics',
             sub_txt='reliability, shortage')

# Arrows
arrow(cx, Y[0]-BH/2, cx, Y[1]+BH/2,
      lbl='LHS design', lbl_color=GOLD)
arrow(cx, Y[1]-BH/2, cx, Y[2]+BH/2)
arrow(cx, Y[2]-BH/2, cx, Y[3]+BH/2)
arrow(cx, Y[3]-BH/2, cx, Y[4]+BH/2)

# Analysis box (below Y[4])
Y_AN_A = Y[4] - BH/2 - 0.45
box_and_text(cx, Y_AN_A, BW, BH-0.08, PURPLE,
             'Scenario discovery\nin parameter space', None, main_fs=6.3)
arrow(cx, Y[4]-BH/2, cx, Y_AN_A+BH/2)

# Annotation at panel bottom
ax.text(cx, PAN_BOT + 0.06, 'Hazard coverage: emergent',
        ha='center', va='bottom', fontsize=5.3, color='#666666', style='italic', zorder=5)

# ── PANEL B — Bottom-up vulnerability ─────────────────────────────────────────
PI = 1
cx = PCX[PI]

box_and_text(cx, Y[0], BW, BH, BLUE,
             'Climate state space /\nForcing variable range',
             sub_txt='ΔP, ΔT, Δdemand grid')
box_and_text(cx, Y[1], BW, BH, ORANGE,
             'Weather generator /\nperturbed scenarios', None)
box_and_text(cx, Y[2], BW, BH+0.04, GRAY,
             'Perturbed scenario ensemble',
             sub_txt='(forcing grid)', text_color=DKGRAY, ec='#909090', lw=0.6)
box_and_text(cx, Y[3], BW, BH, ORANGE,
             'Simulation model', sub_txt='(Pywr-DRB)')

# Outcome box with threshold indicator
box_and_text(cx, Y[4], BW, BH, GREEN,
             'Climate response function',
             sub_txt='performance vs. forcing state')
# Small threshold marker inside green box
ax.plot([cx - BW/2 + 0.2, cx + BW/2 - 0.2],
        [Y[4] - 0.06, Y[4] - 0.06],
        color='white', lw=0.9, linestyle='--', zorder=5, alpha=0.7)
ax.text(cx + BW/2 - 0.25, Y[4] - 0.19,
        '— failure\n   threshold', ha='left', va='center',
        fontsize=4.5, color='white', style='italic', zorder=6)

# Analysis box
Y_AN_B = Y[4] - BH/2 - 0.45
box_and_text(cx, Y_AN_B, BW, BH-0.08, PURPLE,
             'Failure boundary\nin forcing space', None, main_fs=6.3)

# Arrows
arrow(cx, Y[0]-BH/2, cx, Y[1]+BH/2,
      lbl='Grid / stress-test design', lbl_color=GOLD)
arrow(cx, Y[1]-BH/2, cx, Y[2]+BH/2)
arrow(cx, Y[2]-BH/2, cx, Y[3]+BH/2)
arrow(cx, Y[3]-BH/2, cx, Y[4]+BH/2)
arrow(cx, Y[4]-BH/2, cx, Y_AN_B+BH/2)

# Feedback dashed arrow: climate response function → forcing design
# Goes from right side of green box up to the right side of box1, curves right
x_fb = PLX[PI] + PW - 0.08
arrow(x_fb, Y[4], x_fb, Y[0],
      color='#888888', lw=0.9, dashed=True, rad=-0.0,
      lbl='identifies\nfailure region', lbl_side='right', lbl_color='#666666')

ax.text(cx, PAN_BOT + 0.06, 'Hazard coverage: boundary-focused',
        ha='center', va='bottom', fontsize=5.3, color='#666666', style='italic', zorder=5)

# ── PANEL C — Library-and-subsample ──────────────────────────────────────────
PI = 2
cx = PCX[PI]

box_and_text(cx, Y[0], BW, BH, BLUE,
             'Historical streamflow\nrecord', sub_txt='observed data')
box_and_text(cx, Y[1], BW, BH, ORANGE,
             'Kirsch-Nowak generator\n(unsteered)', None)
box_and_text(cx, Y[2], BW, BH+0.04, GRAY,
             'Large library  (10,000+ traces)',
             sub_txt='input space broadly covered', text_color=DKGRAY,
             ec='#909090', lw=0.6)

# Subsampling box (between library and efficient subsample)
Y_SUBSAMP = Y[2] - BH/2 - 0.45
box_and_text(cx, Y_SUBSAMP, BW, BH-0.10, '#7F6000',
             'cLHS subsample\n(input-space coverage)', None,
             main_fs=6.3, text_color=WHITE)

box_and_text(cx, Y[3], BW, BH, ORANGE,
             'Simulation model', sub_txt='(N runs only)')
box_and_text(cx, Y[4], BW, BH, GREEN,
             'Performance metrics', sub_txt='reliability, shortage')

Y_AN_C = Y[4] - BH/2 - 0.45
box_and_text(cx, Y_AN_C, BW, BH-0.08, PURPLE,
             'Robustness ranking /\nscenario discovery', None, main_fs=6.3)

# Arrows
arrow(cx, Y[0]-BH/2, cx, Y[1]+BH/2)
arrow(cx, Y[1]-BH/2, cx, Y[2]+BH/2)
arrow(cx, Y[2]-BH/2, cx, Y_SUBSAMP+(BH-0.10)/2,
      lbl='cLHS design\n(input space)', lbl_color=GOLD)
arrow(cx, Y_SUBSAMP-(BH-0.10)/2, cx, Y[3]+BH/2)
arrow(cx, Y[3]-BH/2, cx, Y[4]+BH/2)
arrow(cx, Y[4]-BH/2, cx, Y_AN_C+BH/2)

ax.text(cx, PAN_BOT + 0.06, 'Hazard coverage: input-space projected',
        ha='center', va='bottom', fontsize=5.3, color='#666666', style='italic', zorder=5)

# ── PANEL D — MOEA-FIND ───────────────────────────────────────────────────────
PI = 3

# Panel D uses a two-column layout:
#   Left column (MOEA box):  PLX[3] + 0.08  to  PLX[3] + 1.25  → width 1.17
#   Right column (main flow): PLX[3] + 1.38  to  PLX[3]+PW-0.10 → width ≈ 1.95
MBW  = 1.17   # MOEA box width
MBH  = 1.65   # MOEA box height (tall — spans generator-to-green vertically)
MBX  = PLX[PI] + 0.08 + MBW/2          # MOEA box center x
MBY  = (Y[1] + Y[3]) / 2              # MOEA box center y (midpoint of gen and green)

MAIN_W  = PW - 0.10 - 1.38 - 0.10     # ≈ 1.95 in
MAIN_CX = PLX[PI] + 1.38 + MAIN_W/2  # main flow center x ≈ PLX[3]+2.36

# Box y for panel D (fewer boxes, spread a bit)
YD = [
    Y[0],                 # YD[0] = Blue (historical data)
    Y[0] - 0.95,          # YD[1] = Generator
    Y[1] - 0.95,          # YD[2] = Green hazard
    Y[2] - 0.95,          # YD[3] = Purple ensemble
]

# Recalculate properly
YD[0] = Y[0]
YD[1] = Y[0] - 0.95
YD[2] = Y[0] - 0.95 - 1.05
YD[3] = Y[0] - 0.95 - 1.05 - 1.00

# Blue input box (full width, at top)
full_cx_D = PLX[PI] + PW/2
box_and_text(full_cx_D, YD[0], BW, BH, BLUE,
             'Historical streamflow\nrecord', sub_txt='observed data')

# Generator box (main column)
box_and_text(MAIN_CX, YD[1], MAIN_W, BH, ORANGE,
             'Kirsch-Nowak\ngenerator', None, main_fs=6.3)

# Green hazard box (main column)
box_and_text(MAIN_CX, YD[2], MAIN_W, BH, GREEN,
             'Drought hazard\ncharacteristics',
             sub_txt='D₁, D₂, D₃', main_fs=6.3)

# Purple ensemble box (main column, full width)
YD_PUR = YD[2] - BH/2 - 0.50
box_and_text(full_cx_D, YD_PUR, BW, BH, PURPLE,
             'Feasible hazard region discovered\n→  Structured ensemble delivered',
             None, main_fs=6.3)

# MOEA box (left column, tall, deep red)
rbox(MBX, MBY, MBW, MBH, RED, ec='#800000', lw=1.2, zorder=3)
txt(MBX, MBY + 0.30, 'Borg MOEA', fs=6.4, fc=WHITE)
txt(MBX, MBY + 0.02, 'ε-dominance\narchive', fs=5.8, fc=WHITE, fw='normal')
txt(MBX, MBY - 0.35, '→  steered\n    decision\n    variables', fs=5.4, fc='#FFCCCC',
    fw='normal', style='italic')

# ── Panel D arrows ────────────────────────────────────────────────────────────

# 1. Blue → Generator (forward, black solid, from center)
arrow(full_cx_D, YD[0]-BH/2, MAIN_CX, YD[1]+BH/2,
      color='black', lw=1.1, rad=-0.12)

# 2. Generator → Green (forward, black solid)
arrow(MAIN_CX, YD[1]-BH/2, MAIN_CX, YD[2]+BH/2, color='black', lw=1.1)

# 3. Green → Purple (forward, black solid, from full-width)
arrow(MAIN_CX, YD[2]-BH/2, full_cx_D, YD_PUR+BH/2,
      color='black', lw=1.1, rad=0.15)

# 4. RED DASHED: Green → MOEA (reversed flow — hazard coverage target)
#    From left side of green box to right side of MOEA box
GN_LEFT_X = MAIN_CX - MAIN_W/2
GN_Y      = YD[2]
MOEA_RIGHT_X = MBX + MBW/2
MOEA_Y_MID   = MBY + 0.20
arrow(GN_LEFT_X, GN_Y, MOEA_RIGHT_X, MOEA_Y_MID,
      color=RED, lw=1.5, dashed=True, rad=0.0,
      lbl='DESIGN\n(desired\ncoverage)', lbl_side='left',
      lbl_color=GOLD)

# 5. RED DASHED: MOEA → Generator (MOEA steers DVs)
#    From top of MOEA box to left side of generator box
GEN_LEFT_X = MAIN_CX - MAIN_W/2
GEN_Y      = YD[1]
MOEA_TOP_X = MBX
MOEA_TOP_Y = MBY + MBH/2
arrow(MOEA_TOP_X, MOEA_TOP_Y, GEN_LEFT_X, GEN_Y,
      color=RED, lw=1.5, dashed=True, rad=-0.3)

# Annotation
ax.text(full_cx_D, PAN_BOT + 0.06, 'Hazard coverage: structured by construction',
        ha='center', va='bottom', fontsize=5.3, color='#006020',
        fontweight='bold', style='italic', zorder=5)

# ── Shared hazard-space coverage strip ────────────────────────────────────────
rng = np.random.default_rng(42)

strip_labels = ['uncontrolled', 'boundary-focused', 'input-projected', 'structured']
strip_colors = ['#888888', '#606060', '#606060', '#1B5E20']

for pi in range(4):
    sl = PLX[pi]
    sr = sl + PW
    # Strip background
    ax.add_patch(plt.Rectangle((sl, STRIP_BOT), PW, STRIP_H,
                                fc=SLGN, ec='#B0C8A8', lw=0.6, zorder=2))
    # x and y axis labels (tiny)
    ax.text(sl + PW/2, STRIP_BOT + STRIP_H - 0.11,
            'Drought hazard space  (severity → , duration ↑)',
            ha='center', va='top', fontsize=5.0, color='#336622',
            style='italic', zorder=5)

    # Dot patterns
    sx_range = (sl + 0.18, sl + PW - 0.18)
    sy_range = (STRIP_BOT + 0.10, STRIP_BOT + STRIP_H - 0.20)
    sw = sx_range[1] - sx_range[0]
    sh = sy_range[1] - sy_range[0]

    if pi == 0:   # Panel A: random, clustered toward center
        n = 32
        x_raw = rng.normal(0.5, 0.22, n)
        y_raw = rng.normal(0.5, 0.22, n)
        x_raw = np.clip(x_raw, 0.05, 0.95)
        y_raw = np.clip(y_raw, 0.05, 0.95)
        xs = sx_range[0] + x_raw * sw
        ys = sy_range[0] + y_raw * sh
        sizes = rng.uniform(6, 14, n)
        ax.scatter(xs, ys, s=sizes, c='#4472C4', alpha=0.65, zorder=4, linewidths=0)

    elif pi == 1:  # Panel B: concentrated near a diagonal threshold line
        n = 28
        # Dots near diagonal y = x (the failure boundary)
        t = rng.uniform(0.1, 0.9, n)
        perp = rng.normal(0, 0.06, n)
        x_raw = np.clip(t + perp, 0.05, 0.95)
        y_raw = np.clip(t - perp, 0.05, 0.95)
        xs = sx_range[0] + x_raw * sw
        ys = sy_range[0] + y_raw * sh
        # Color by failure side
        failure = (x_raw + y_raw) > 1.0
        colors  = ['#C00000' if f else '#4472C4' for f in failure]
        ax.scatter(xs, ys, s=10, c=colors, alpha=0.70, zorder=4, linewidths=0)
        # Draw the threshold boundary
        bx = np.array([sx_range[0], sx_range[1]])
        by = np.array([sy_range[0] + sh, sy_range[0]])  # diagonal from top-left to bottom-right
        ax.plot(bx, by, color='#333333', lw=0.9, linestyle='--', zorder=4)
        ax.text(sx_range[0] + 0.18, sy_range[0] + sh * 0.85,
                'failure\nboundary', ha='center', va='center',
                fontsize=4.3, color='#333333', style='italic', zorder=5)

    elif pi == 2:  # Panel C: moderate structure, nonlinear projection distortion
        n = 30
        # Slightly structured but with nonlinear bunching toward low-severity
        x_raw = rng.beta(1.5, 2.5, n)   # skewed toward low severity
        y_raw = rng.beta(1.5, 2.0, n)
        xs = sx_range[0] + x_raw * sw
        ys = sy_range[0] + y_raw * sh
        ax.scatter(xs, ys, s=9, c='#6B6B6B', alpha=0.65, zorder=4, linewidths=0)

    else:          # Panel D: near-uniform grid (structured MOEA archive)
        n_side = 6
        gx = np.linspace(0.08, 0.92, n_side)
        gy = np.linspace(0.10, 0.90, n_side)
        # Restrict to roughly elliptical feasible region
        pts = []
        for gxi in gx:
            for gyi in gy:
                # Elliptical constraint (feasible region is not a full square)
                if ((gxi - 0.5)/0.45)**2 + ((gyi - 0.5)/0.45)**2 <= 1.05:
                    pts.append((gxi, gyi))
        pts = np.array(pts)
        # Add tiny jitter for realism
        pts += rng.normal(0, 0.015, pts.shape)
        xs = sx_range[0] + pts[:, 0] * sw
        ys = sy_range[0] + pts[:, 1] * sh
        ax.scatter(xs, ys, s=11, c='#1B7A3A', alpha=0.80, zorder=4,
                   marker='D', linewidths=0)
        # Draw feasible region boundary (ellipse)
        theta = np.linspace(0, 2*np.pi, 80)
        ell_x = sx_range[0] + (0.5 + 0.45*np.cos(theta)) * sw
        ell_y = sy_range[0] + (0.5 + 0.45*np.sin(theta)) * sh
        ax.plot(ell_x, ell_y, color='#1B5E20', lw=0.8, linestyle=':', zorder=4, alpha=0.6)

    # Strip annotation (bottom right corner of each strip cell)
    ax.text(sl + PW - 0.08, STRIP_BOT + 0.06,
            strip_labels[pi], ha='right', va='bottom',
            fontsize=5.0, color=strip_colors[pi], fontweight='bold', zorder=5)

# ── Vertical dividers between panels in the strip ────────────────────────────
for pi in range(1, 4):
    xd = PLX[pi] - PGAP/2
    ax.plot([xd, xd], [STRIP_BOT, STRIP_BOT + STRIP_H],
            color='#C8D8C0', lw=0.6, zorder=3)

# ── Strip title (left side) ───────────────────────────────────────────────────
ax.text(LM - 0.04, STRIP_BOT + STRIP_H/2,
        'Coverage\nin drought\nhazard space',
        ha='right', va='center', fontsize=5.5,
        color='#336622', fontweight='bold',
        multialignment='center', zorder=5)

# ── Legend for strip dot colors ───────────────────────────────────────────────
leg_x = PLX[1] + 0.10
leg_y = STRIP_BOT + 0.06
ax.scatter([leg_x + 0.10], [leg_y + 0.12], s=10, c='#4472C4', alpha=0.75,
           zorder=6, linewidths=0)
ax.text(leg_x + 0.22, leg_y + 0.12, 'pass', ha='left', va='center',
        fontsize=4.5, color='#4472C4', zorder=6)
ax.scatter([leg_x + 0.55], [leg_y + 0.12], s=10, c='#C00000', alpha=0.75,
           zorder=6, linewidths=0)
ax.text(leg_x + 0.67, leg_y + 0.12, 'fail', ha='left', va='center',
        fontsize=4.5, color='#C00000', zorder=6)
ax.text(leg_x, leg_y + 0.12, '(panel b only:)', ha='left', va='center',
        fontsize=4.5, color='#555555', zorder=6)

# ── DESIGN legend box (top-right corner, outside panels) ─────────────────────
LEGEND_X = PLX[3] + PW + 0.03
LEGEND_Y = PAN_BOT + 2.2
if LEGEND_X + 0.12 < FW:
    pass   # no room outside panels — skip external legend

# ── Arrow legend inside panel A (small, bottom area) ─────────────────────────
leg_ax_x = PLX[0] + 0.12
leg_ax_y = PAN_BOT + 0.62
ax.annotate('', xy=(leg_ax_x + 0.38, leg_ax_y), xytext=(leg_ax_x, leg_ax_y),
            arrowprops=dict(arrowstyle='->', color='black', lw=0.9), zorder=6)
ax.text(leg_ax_x + 0.42, leg_ax_y, 'forward flow', ha='left', va='center',
        fontsize=4.5, color='black', zorder=6)

leg_rx_x = PLX[0] + 0.12
leg_rx_y = PAN_BOT + 0.40
ax.annotate('', xy=(leg_rx_x + 0.38, leg_rx_y), xytext=(leg_rx_x, leg_rx_y),
            arrowprops=dict(arrowstyle='->', color=RED, lw=0.9,
                            linestyle=(0, (4, 2.5))), zorder=6)
ax.text(leg_rx_x + 0.42, leg_rx_y, 'search direction (D only)',
        ha='left', va='center', fontsize=4.5, color=RED, zorder=6)

# Badge legend
bx = leg_ax_x
by = PAN_BOT + 0.20
ax.text(bx + 0.12, by, 'DESIGN',
        ha='center', va='center', fontsize=4.5, color=GOLD, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.10', fc=WHITE, ec=GOLD, lw=0.7),
        zorder=6)
ax.text(bx + 0.50, by, '= structured design step',
        ha='left', va='center', fontsize=4.5, color='#555555', zorder=6)

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir = Path(__file__).parent
fig.savefig(out_dir / 'figure1_mockup.svg', format='svg',
            bbox_inches='tight', dpi=150)
fig.savefig(out_dir / 'figure1_mockup.pdf', format='pdf',
            bbox_inches='tight')
print('Saved: figure1_mockup.svg  +  figure1_mockup.pdf')
plt.close(fig)
