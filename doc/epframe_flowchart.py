#!/usr/bin/env python3
"""
epframe_flowchart.py
Two-page PDF:
  Page 1 — Solution flowchart
  Page 2 — Variable definitions table

Uses matplotlib mathtext ($...$) for all mathematical symbols —
no LaTeX installation required, no font encoding issues.

Layout constants are collected near the top of each page function
so they are easy to adjust.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

# ── shared colours ─────────────────────────────────────────────────────────────
C_INIT  = '#D6E4F0'   # light blue   – initialisation / output
C_LOOP  = '#D5E8D4'   # light green  – main loop boxes
C_SUB   = '#FFF2CC'   # yellow       – active-set sub-loop
C_TERM  = '#F8CECC'   # light red    – termination
C_HEAD  = '#1F3864'   # dark blue    – header / start / end / borders
C_WARN  = '#D79B00'   # amber        – sub-loop border
C_WHITE = 'white'
C_BLACK = 'black'

# ── helper: draw a box with wrapped centred text ──────────────────────────────
def box(ax, x, y, w, h, facecolor, lines, fontsize=8.0,
        bold_first=False, edgecolor=C_HEAD, lw=0.9, radius=0.008, zorder=3):
    """
    Draw a rounded rectangle at (x, y) with width w, height h (all in axes
    fraction units), containing multi-line centred text.
    (x, y) is the BOTTOM-LEFT corner.
    """
    patch = FancyBboxPatch((x, y), w, h,
                           boxstyle=f"round,pad=0,rounding_size={radius}",
                           facecolor=facecolor, edgecolor=edgecolor,
                           linewidth=lw, zorder=zorder,
                           transform=ax.transAxes, clip_on=False)
    ax.add_patch(patch)
    n = len(lines)
    for i, line in enumerate(lines):
        ty = y + h/2 + (n/2 - i - 0.5) * (fontsize + 1.5) / (ax.get_figure().get_figheight() * 72)
        weight = 'bold' if (bold_first and i == 0) else 'normal'
        ax.text(x + w/2, ty, line,
                ha='center', va='center', fontsize=fontsize,
                fontweight=weight, color=C_BLACK,
                transform=ax.transAxes, clip_on=False, zorder=zorder+1)

def arrow_v(ax, x, y_top, length, color=C_HEAD, lw=1.0):
    """Downward arrow in axes-fraction coords."""
    ax.annotate('', xy=(x, y_top - length), xytext=(x, y_top),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color=color, lw=lw),
                zorder=5)

def bracket_loop(ax, x_right, y_bot, y_top, x_re_entry, color=C_HEAD):
    """Right-side loop-back bracket from y_bot back up to y_top."""
    # vertical line on right
    ax.plot([x_right, x_right], [y_bot, y_top],
            transform=ax.transAxes, color=color, lw=1.0, ls='--', zorder=2)
    # horizontal line at top back to entry
    ax.annotate('', xy=(x_re_entry, y_top), xytext=(x_right, y_top),
                xycoords='axes fraction', textcoords='axes fraction',
                arrowprops=dict(arrowstyle='->', color=color, lw=1.0,
                                linestyle='dashed'),
                zorder=5)
    # horizontal at bottom
    ax.plot([x_right - 0.01, x_right], [y_bot, y_bot],
            transform=ax.transAxes, color=color, lw=1.0, ls='--', zorder=2)

# ══════════════════════════════════════════════════════════════════════════════
# PAGE 1 — FLOWCHART
# ══════════════════════════════════════════════════════════════════════════════
def page_flowchart(pdf):
    fig = plt.figure(figsize=(8.5, 11))
    ax  = fig.add_axes([0, 0, 1, 1])   # full-page axes
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')

    # ── layout constants (edit these to reposition / resize) ─────────────────
    LM   = 0.25   # left margin (axes fraction)
    RM   = 0.75   # right margin
    CX   = 0.50   # horizontal centre

    BW_FULL = RM - LM          # full-width box width  ≈ 0.86
    BW_LOOP = BW_FULL * 0.85   # inner loop box width  (−15 %)
    LM_LOOP = LM + (BW_FULL - BW_LOOP) / 2   # re-centred left edge

    BH_SM  = 0.032  # small  box height
    BH_STD = 0.038  # standard box height  (×1.15 vs original)
    BH_LG  = 0.055  # large  box height    (×1.15 vs original)
    BH_XL  = 0.075  # extra-large (active-set sub-loop)
    BH_TERM= 0.090  # termination block

    GAP  = 0.012   # vertical gap between boxes
    ARR  = 0.020   # arrow length

    FS_BOX = 7.8   # font size inside boxes
    FS_SUB = 7.4   # font size in sub-loop

    # ── title ─────────────────────────────────────────────────────────────────
    ax.text(CX, 0.975,
            'Incremental Elastic-Plastic Frame Analysis — Solution Algorithm',
            ha='center', va='top', fontsize=11, fontweight='bold',
            color=C_HEAD, transform=ax.transAxes)

    y = 0.945   # running y: TOP of next box to draw

    # ── START ─────────────────────────────────────────────────────────────────
    sw = 0.10
    bx_s = CX - sw/2
    box(ax, bx_s, y - BH_SM, sw, BH_SM, C_HEAD, ['START'],
        fontsize=10, bold_first=True, edgecolor=C_HEAD)
    ax.texts[-1].set_color(C_WHITE)
    y -= BH_SM + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── INITIALISATION ────────────────────────────────────────────────────────
    init_lines = [
        'Build $\\mathbf{C}$  (compatibility matrix from frame geometry)',
        'Initialise $\\mathbf{SF}$, $\\mathbf{SA}$ ;   '
        '$\\mathbf{CM} = \\mathbf{CT} = \\mathbf{CD} = 0$',
        'Compute  $DI = 3 \\cdot NE - ND$  and  '
        '$DLMT = 0.1 \\times \\max(\\mathrm{frame\\ dim.})$',
    ]
    box(ax, LM, y - BH_LG, BW_FULL, BH_LG, C_INIT, init_lines,
        fontsize=FS_BOX, bold_first=True)
    y -= BH_LG + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── outer loop dashed border (drawn first, behind everything) ─────────────
    loop_top = y

    # ── stiffness assembly ────────────────────────────────────────────────────
    stiff_lines = [
        'Form $\\mathbf{S}$ from current $\\mathbf{SF}$ (bending) and $\\mathbf{SA}$ (axial)',
        '$\\mathbf{K}_e = \\mathbf{C}^T \\cdot \\mathbf{S} \\cdot \\mathbf{C}$'
        '          (elastic stiffness)',
        '$\\mathbf{K}_g = $ assemble\\_geometric($\\mathbf{CT}, \\mathbf{C}, L$)'
        '   (geometric stiffness)',
        '$\\mathbf{K}_{tot} = \\mathbf{K}_e + \\mathbf{K}_g$',
    ]
    box(ax, LM_LOOP, y - BH_LG, BW_LOOP, BH_LG, C_LOOP, stiff_lines,
        fontsize=FS_BOX, bold_first=True)
    y -= BH_LG + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── active-set sub-loop ───────────────────────────────────────────────────
    sub_lines = [
        'ACTIVE-SET SOLVE',
        'Solve  $\\mathbf{K}_{tot} \\cdot \\delta = \\mathbf{P}$  for free DOFs',
        'Lagrange multipliers $\\Rightarrow$ reaction forces at constrained DOFs',
        'If any unidirectional reaction has wrong sign:',
        '  release that constraint $\\rightarrow$ update active set',
        'Repeat until all reactions are sign-consistent',
    ]
    sub_lm = LM_LOOP + 0.015
    sub_bw = BW_LOOP - 0.03
    box(ax, sub_lm, y - BH_XL, sub_bw, BH_XL, C_SUB, sub_lines,
        fontsize=FS_SUB, bold_first=True, edgecolor=C_WARN, lw=1.3)
    # dashed amber border overlay
    sub_patch = FancyBboxPatch((sub_lm, y - BH_XL), sub_bw, BH_XL,
                               boxstyle='round,pad=0,rounding_size=0.006',
                               facecolor='none', edgecolor=C_WARN,
                               linewidth=1.5, linestyle='--', zorder=6,
                               transform=ax.transAxes, clip_on=False)
    ax.add_patch(sub_patch)
    y -= BH_XL + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── forces + checks ───────────────────────────────────────────────────────
    fc_lines = [
        '$\\mathbf{F} = \\mathbf{S} \\cdot \\mathbf{C}^T \\cdot \\delta$'
        '                (member force rates: moments + axial)',
        'Check $|CT[el]| \\geq Py[el]$  $\\Rightarrow$  COMPRESSION YIELD  (stop)',
        'If $NCYCL > DI$              $\\Rightarrow$  MECHANISM WARNING',
    ]
    box(ax, LM_LOOP, y - BH_LG, BW_LOOP, BH_LG, C_LOOP, fc_lines,
        fontsize=FS_BOX)
    y -= BH_LG + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── P-M quadratic ─────────────────────────────────────────────────────────
    pm_lines = [
        'Solve for load factor $\\alpha$  (P–M quadratic at each element end)',
        '$A\\alpha^2 + B\\alpha + C_{coef} = 0$'
        '      take smallest positive real root',
        '$A = M_p \\dot{p}^2 / P_y^2$'
        '          $B = s\\,\\dot{m} + 2 M_p P_0 \\dot{p} / P_y^2$',
        '$C_{coef} = s\\,M_0 - M_p\\,(1 - (P_0/P_y)^2)$',
    ]
    box(ax, LM_LOOP, y - BH_LG, BW_LOOP, BH_LG, C_LOOP, pm_lines,
        fontsize=FS_BOX)
    y -= BH_LG + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── cumulative update ─────────────────────────────────────────────────────
    upd_lines = [
        'Update:  $\\mathbf{CD} {+}= \\alpha\\delta$,  '
        '$\\mathbf{CM} {+}= \\alpha\\dot{m}$,  '
        '$\\mathbf{CT} {+}= \\alpha\\dot{p}$,  '
        '$\\lambda {+}= \\alpha$',
        'Check  $\\max|\\mathbf{CD}| > DLMT$  $\\Rightarrow$  DISPLACEMENT LIMIT  (stop)',
    ]
    box(ax, LM_LOOP, y - BH_STD, BW_LOOP, BH_STD, C_LOOP, upd_lines,
        fontsize=FS_BOX)
    y -= BH_STD + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── hinge modification ────────────────────────────────────────────────────
    hinge_lines = [
        'Modify $\\mathbf{SF}$ at new hinge:',
        'Near-end stiffness $\\rightarrow 0$  (moment released)',
        'Far-end  stiffness $\\rightarrow 0.75 \\times$ prior value'
        '   ($4EI/L \\rightarrow 3EI/L$)',
    ]
    box(ax, LM_LOOP, y - BH_LG, BW_LOOP, BH_LG, C_LOOP, hinge_lines,
        fontsize=FS_BOX, bold_first=True)
    y -= BH_LG

    loop_bot = y   # bottom of loop region

    # draw outer loop dashed rectangle
    loop_rect = FancyBboxPatch((LM, loop_bot - 0.008),
                               BW_FULL, loop_top - loop_bot + 0.008,
                               boxstyle='round,pad=0,rounding_size=0.008',
                               facecolor='none', edgecolor=C_HEAD,
                               linewidth=1.2, linestyle=(0,(6,3)),
                               zorder=2, transform=ax.transAxes, clip_on=False)
    ax.add_patch(loop_rect)
    ax.text(LM - 0.005, (loop_top + loop_bot)/2,
            'LOAD INCREMENT LOOP',
            ha='center', va='center', fontsize=7, fontweight='bold',
            color=C_HEAD, rotation=90, transform=ax.transAxes)

    # loop-back arrow on right side
    bracket_loop(ax, RM - 0.02, loop_bot + 0.02, loop_top - 0.02,
                 CX + BW_LOOP/2, color=C_HEAD)

    y = loop_bot - 0.008 - ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── termination ───────────────────────────────────────────────────────────
    term_lines = [
        'TERMINATION  (first condition met stops the analysis)',
        '$\\mathbf{K}_{tot}$ singular or negative eigenvalue'
        '  $\\Rightarrow$  buckling / mechanism',
        '$|CT| \\geq Py$ in any element'
        '                $\\Rightarrow$  compression yield',
        '$\\max|\\mathbf{CD}| > DLMT$'
        '                      $\\Rightarrow$  displacement limit exceeded',
        '$NCYCL > DI$'
        '                             $\\Rightarrow$  kinematic mechanism likely',
        '$\\alpha > 10^9$'
        '                               $\\Rightarrow$  no further hinge possible',
    ]
    box(ax, LM, y - BH_TERM, BW_FULL, BH_TERM, C_TERM, term_lines,
        fontsize=FS_BOX, bold_first=True)
    y -= BH_TERM + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── output ────────────────────────────────────────────────────────────────
    out_lines = [
        'Write output files  (.out  and  .csv)',
        'Generate plots  (epframe\\_viz.py)',
    ]
    box(ax, LM, y - BH_STD, BW_FULL, BH_STD, C_INIT, out_lines,
        fontsize=FS_BOX)
    y -= BH_STD + ARR
    arrow_v(ax, CX, y + ARR, ARR)

    # ── END ───────────────────────────────────────────────────────────────────
    ew = 0.10
    box(ax, CX - ew/2, y - BH_SM, ew, BH_SM, C_HEAD, ['END'],
        fontsize=10, bold_first=True, edgecolor=C_HEAD)
    ax.texts[-1].set_color(C_WHITE)

    # footer
    ax.text(CX, 0.012,
            'epframe.py — H. P. Gavin, Duke University CEE  |  '
            'Incremental Elastic-Plastic Frame Analysis',
            ha='center', va='bottom', fontsize=6.5,
            color='grey', transform=ax.transAxes, style='italic')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# PAGE 2 — VARIABLE DEFINITIONS TABLE
# ══════════════════════════════════════════════════════════════════════════════
def page_variables(pdf):
    fig = plt.figure(figsize=(8.5, 11))
    ax  = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')

    # ── layout constants ──────────────────────────────────────────────────────
    LM      = 0.06    # left margin
    RM      = 0.94    # right margin
    TW      = RM - LM # total table width
    SYM_W   = 0.16    # symbol column width  (fraction of total width)
    DEF_W   = 1.0 - SYM_W   # definition column width (reduced 20 %)
    ROW_H   = 0.024   # row height
    FS      = 8.0     # font size
    HEAD_Y  = 0.945   # top of table

    # ── title ─────────────────────────────────────────────────────────────────
    ax.text(0.50, 0.975, 'Variable Definitions — Alphabetical',
            ha='center', va='top', fontsize=11, fontweight='bold',
            color=C_HEAD, transform=ax.transAxes)

    # ── table rows: (symbol_mathtext, definition_text) ───────────────────────
    rows = [
        # header
        ('Symbol', 'Definition'),
        # data
        (r'$\alpha$',
         r'Load factor increment to the next plastic hinge (scalar solved each step)'),
        (r'$A,\, B,\, C_{coef}$',
         r'Coefficients of the quadratic equation for $\alpha$ in P-M hinge detection'),
        (r'$\mathbf{C}$',
         r'Compatibility matrix — maps member deformations to structural DOF displacements'),
        (r'$\mathbf{CD}$',
         r'Cumulative nodal displacement vector (updated each increment)'),
        (r'$\mathbf{CM}$',
         r'Cumulative member end-moment vector (updated each increment)'),
        (r'$\mathbf{CT}$',
         r'Cumulative member axial tension force vector (updated each increment)'),
        (r'$\delta$',
         r'Nodal displacement rate vector — solution of $\mathbf{K}_{tot}\cdot\delta = \mathbf{P}$'),
        (r'$DI$',
         r'Degree of static indeterminacy $= 3\cdot NE - ND$'),
        (r'$DLMT$',
         r'Displacement limit $= 0.1 \times \max(\mathrm{frame\ dimension})$'),
        (r'$E$',
         r'Elastic modulus'),
        (r'$\mathbf{F}$',
         r'Member force rate vector (moments + axial) $= \mathbf{S}\cdot\mathbf{C}^T\cdot\delta$'),
        (r'$F_y$',
         r'Yield stress'),
        (r'$\mathbf{K}_e$',
         r'Global elastic stiffness matrix $= \mathbf{C}^T \cdot \mathbf{S} \cdot \mathbf{C}$'),
        (r'$\mathbf{K}_g$',
         r'Global geometric stiffness matrix (element $\tilde{G}$ matrices weighted by $\mathbf{CT}$)'),
        (r'$\mathbf{K}_{tot}$',
         r'Total stiffness matrix $= \mathbf{K}_e + \mathbf{K}_g$'),
        (r'$\lambda$',
         r'Cumulative load factor (sum of all $\alpha$ increments)'),
        (r'$L$',
         r'Element length'),
        (r'$M_0$',
         r'Cumulative moment at an element end at the start of the current increment'),
        (r'$\dot{m}$',
         r'Moment rate at an element end (from $\mathbf{F}$)'),
        (r'$M_p$',
         r'Plastic moment capacity $= Z \cdot F_y$'),
        (r'$M_u$',
         r'Reduced plastic moment under axial load $= M_p\,(1-(P/P_y)^2)$'),
        (r'$ND$',
         r'Number of free (unconstrained) structural degrees of freedom'),
        (r'$NE$',
         r'Number of frame elements'),
        (r'$\mathbf{P}$',
         r'Applied nodal load vector'),
        (r'$P_0$',
         r'Cumulative axial force in an element at the start of the current increment'),
        (r'$\dot{p}$',
         r'Axial force rate in an element (from $\mathbf{F}$)'),
        (r'$P_y$',
         r'Axial yield force $= A \cdot F_y$'),
        (r'$s$',
         r'Sign of $M_0$ (or of $\dot{m}$ when $M_0\approx 0$) — side of interaction surface'),
        (r'$\mathbf{S}$',
         r'Block-diagonal member stiffness matrix ($\mathbf{SF}$ bending + $\mathbf{SA}$ axial blocks)'),
        (r'$\mathbf{SA}$',
         r'Axial stiffness coefficients $= EA/L$ per element'),
        (r'$\mathbf{SF}$',
         r'Bending stiffness coefficients ($4EI/L$, $2EI/L$ terms; reduced at hinge locations)'),
        (r'$Z$',
         r'Plastic section modulus'),
    ]

    # ── draw header ───────────────────────────────────────────────────────────
    y = HEAD_Y
    # header background
    hdr = FancyBboxPatch((LM, y - ROW_H), TW, ROW_H,
                         boxstyle='square,pad=0',
                         facecolor=C_HEAD, edgecolor=C_HEAD,
                         zorder=3, transform=ax.transAxes, clip_on=False)
    ax.add_patch(hdr)
    ax.text(LM + TW * SYM_W / 2, y - ROW_H/2, rows[0][0],
            ha='center', va='center', fontsize=FS, fontweight='bold',
            color='white', transform=ax.transAxes, zorder=4)
    ax.text(LM + TW * SYM_W + TW * DEF_W / 2, y - ROW_H/2, rows[0][1],
            ha='center', va='center', fontsize=FS, fontweight='bold',
            color='white', transform=ax.transAxes, zorder=4)
    y -= ROW_H

    # ── draw data rows ────────────────────────────────────────────────────────
    for i, (sym, defn) in enumerate(rows[1:]):
        bg = C_INIT if i % 2 == 0 else 'white'
        row_patch = FancyBboxPatch((LM, y - ROW_H), TW, ROW_H,
                                   boxstyle='square,pad=0',
                                   facecolor=bg, edgecolor='#AAAAAA',
                                   linewidth=0.4, zorder=3,
                                   transform=ax.transAxes, clip_on=False)
        ax.add_patch(row_patch)

        # symbol (mathtext, centred in symbol column)
        ax.text(LM + TW * SYM_W / 2, y - ROW_H/2, sym,
                ha='center', va='center', fontsize=FS, color=C_BLACK,
                transform=ax.transAxes, zorder=4)

        # definition (mathtext, left-aligned in definition column)
        ax.text(LM + TW * SYM_W + 0.008, y - ROW_H/2, defn,
                ha='left', va='center', fontsize=FS - 0.3, color=C_BLACK,
                transform=ax.transAxes, zorder=4)

        y -= ROW_H

    # vertical divider between columns
    ax.plot([LM + TW * SYM_W, LM + TW * SYM_W],
            [y, HEAD_Y],
            transform=ax.transAxes, color='#888888', lw=0.5, zorder=5)

    # outer border
    outer = FancyBboxPatch((LM, y), TW, HEAD_Y - y,
                           boxstyle='square,pad=0',
                           facecolor='none', edgecolor=C_HEAD,
                           linewidth=1.0, zorder=6,
                           transform=ax.transAxes, clip_on=False)
    ax.add_patch(outer)

    # footer
    ax.text(0.50, 0.012,
            'epframe.py — H. P. Gavin, Duke University CEE',
            ha='center', va='bottom', fontsize=6.5,
            color='grey', transform=ax.transAxes, style='italic')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)


# ── MAIN ──────────────────────────────────────────────────────────────────────
out_path = './epframe_flowchart.pdf'
with PdfPages(out_path) as pdf:
    page_flowchart(pdf)
    page_variables(pdf)

print(f'Written: {out_path}')
