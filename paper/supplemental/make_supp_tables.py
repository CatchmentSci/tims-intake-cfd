"""Typeset the Supplemental Materials tables from the analysis CSVs.

Reads cfd/results/supplemental/S1_dimensions.csv, S2_bore_reynolds.csv, S3_meshes.csv and
S4_runs.csv (written by cfd/geometry/supplement_tables.py) and writes
tab_S1.tex ... tab_S4.tex next to this script. Change numbers in the CSVs (via the analysis
script), not here; change only the layout here.  Run: python -I make_supp_tables.py
"""
import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, '..', '..', 'cfd', 'results', 'supplemental'))
VELOCITIES = ['0.154', '0.3', '0.585', '1', '1.25', '1.5']
RAG = r'>{\raggedright\arraybackslash}'


def read(name):
    with open(os.path.join(DATA, name), encoding='utf-8', newline='') as f:
        return list(csv.DictReader(f))


def tex(s):
    """Escape plain CSV text for LaTeX and tidy a few typographic details."""
    s = (s or '').strip()
    s = re.sub(r'~(?=\d)', '\x01', s)
    s = re.sub(r'\^(\d)', '\x02' + r'\1' + '\x03', s)
    s = s.replace('\\', r'\textbackslash{}')
    for a, b in [('&', r'\&'), ('%', r'\%'), ('_', r'\_'), ('#', r'\#'), ('~', r'\textasciitilde{}'),
                 ('^', r'\textasciicircum{}')]:
        s = s.replace(a, b)
    s = re.sub(r' x (?=\d)', r'~$\\times$~', s)
    s = re.sub(r'(?<=\d) - (?=\d)', r' $-$ ', s)
    s = re.sub(r'(?<=\d)x(?=\d)', r'$\\times$', s)
    s = s.replace('\x01', r'$\approx$').replace('\x02', r'$^{').replace('\x03', '}$')
    s = s.replace(' deg', r'$^\circ$').replace('+-', r'$\pm$').replace('>=', r'$\geq$')
    return s


def tex_note(s):
    """tex() for the prose of the Table S3 notes: number ranges and Hagen--Poiseuille with en dashes, and
    a bare R in italics, as in the paper. Not used for Table S1, whose part numbers (D1785-21) need hyphens."""
    s = tex(s)
    s = re.sub(r'(?<=\d)-(?=\d)', '--', s)
    s = s.replace('Hagen-Poiseuille', 'Hagen--Poiseuille')
    s = re.sub(r'(?<![\w$\\])R(?=[ =,.])', '$R$', s)
    return s


STATUS = {'published': 'pub.', 'catalogue': r'cat.$^\dagger$', 'assumed': r'ass.$^\ast$',
          'derived': 'der.'}


# Table S1 source text -> citations (refs 1-2). Applied after
# tex(), so keys match the escaped text. Sources that could not be verified stay as plain text.
S1_CITES = [
    ('Phillips et al. 2000', r'\citet{Phillips2000}'),
    ('Russell et al. 2000', r'\citet{Russell2000}'),
    (r'Smith \& Owens 2014', r'\citet{SmithOwens2014}'),
    ('Elliott et al. 2017', r'\citet{Elliott2017}'),
    ('McDonald et al. 2010', r'\citet{McDonald2010}'),
    ('Foets et al. 2020', r'\citet{Foets2020}'),
    ('ASTM D1785', r'\citet{ASTMD1785}'),
    ('ASTM D2466', r'\citet{ASTMD2466}'),
    ('BS EN 10226', r'\citet{BSEN10226}'),
    ('UF25B datasheet', r'\citet{Cynergy2020}'),
    (r'Nalgene 980, 5/16~$\times$~1/2 in (Fisher 10448612)', r'Nalgene 980, 5/16~$\times$~1/2 in \citep{Nalgene980}'),
    ('Nalgene 4256-0414', r'Nalgene 4256-0414 \citep{NalgeneFunnel}'),
]


def cite_sources(s):
    for a, b in S1_CITES:
        s = s.replace(a, b)
    return s


def table_s1():
    rows = read('S1_dimensions.csv')
    out = [r'\begin{longtable}{@{}' + RAG + r'p{0.13\textwidth}' + RAG + r'p{0.18\textwidth}' + RAG +
           r'p{0.08\textwidth}' + RAG + r'p{0.36\textwidth}' + RAG + r'p{0.12\textwidth}@{}}',
           r"\caption{Dimensions of the modelled designs (mm). Status: published (pub.), builders'"
           r' parts list or manufacturer data (cat.$^\dagger$), assumed (ass.$^\ast$) or derived'
           r' from these (der.). Mounting hardware (uprights, booms, clamps, ice deflectors) is'
           r' omitted for all designs.}\label{tab:S1}\\',
           r'\toprule', r'\textbf{Item} & \textbf{Value} & \textbf{Status} & '
           r'\textbf{Source and reasoning} & \textbf{Range} \\', r'\midrule', r'\endfirsthead',
           r'\toprule', r'\textbf{Item} & \textbf{Value} & \textbf{Status} & '
           r'\textbf{Source and reasoning} & \textbf{Range} \\', r'\midrule', r'\endhead',
           r'\bottomrule', r'\endlastfoot']
    design = None
    for r in rows:
        if r['design'] != design:
            design = r['design']
            if len(out) > 14:
                out.append(r'\addlinespace')
            out.append(r'\multicolumn{5}{@{}l}{\textbf{' + tex(design) + r'}} \\')
        src = '; '.join(x for x in (cite_sources(tex(r['source'])), tex(r['reasoning'])) if x)
        out.append(' & '.join([tex(r['item']), tex(r['value_mm']),
                               r'/\allowbreak{}'.join(STATUS.get(x.strip(), tex(x)) for x in r['status'].split('/')), src,
                               tex(r['range'])]) + r' \\')
    out.append(r'\end{longtable}')
    return '\n'.join(out) + '\n'


def table_s2():
    rows = read('S2_bore_reynolds.csv')
    out = [r'\begin{table}[!htb]', r'\centering',
           r'\caption{Bore Reynolds number $Re_b = 4q/(\pi d \nu)$ of each bore the sampled flow'
           r' passes through, from the medium-mesh sweep ($q = R\,U\,A_\mathrm{in}$; water at'
           r' \SI{20}{\celsius}, $\nu = \SI{1.004e-6}{m^2.s^{-1}}$). Values above 2000 in bold; $Re_b$'
           r' alone does not decide transition in these short, disturbed passages.}\label{tab:S2}',
           r'\small', r'\begin{tabular}{@{}llrrrrrrr@{}}', r'\toprule',
           r'& & & \multicolumn{6}{c}{$U$ (\si{m.s^{-1}})} \\ \cmidrule(l){4-9}',
           r'\textbf{Design} & \textbf{Bore} & $d$ (mm) & ' +
           ' & '.join(VELOCITIES) + r' \\', r'\midrule']
    design = None
    for r in rows:
        name = tex(r['design']) if r['design'] != design else ''
        design = r['design']
        vals = [(r'\textbf{' + r['Re_U' + v] + '}') if float(r['Re_U' + v]) > 2000 else r['Re_U' + v]
                for v in VELOCITIES]
        out.append(' & '.join([name, tex(r['bore']), r['d_mm']] + vals) + r' \\')
    out += [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    return '\n'.join(out) + '\n'


def table_s3():
    rows = read('S3_meshes.csv')
    out = [r'\begin{table}[!htb]', r'\centering',
           r'\caption{Mesh convergence at \SI{0.3}{m.s^{-1}}: scale $s$, cells, faces across each'
           r' bore and $R$ on each mesh, with the effective refinement ratio $r$ from the next coarser mesh,'
           r' $(N_\mathrm{finer}/N_\mathrm{coarser})^{1/3}$, the medium-mesh error, the GCI of the fine mesh'
           r' \citep{Celik2008} and the best estimate $R_\infty$.}\label{tab:S3}',
           r'\footnotesize', r'\setlength{\tabcolsep}{4pt}',
           r'\begin{tabular}{@{}llrrr' + RAG + r'p{0.24\textwidth}rrrr@{}}', r'\toprule',
           r'\textbf{Design} & \textbf{Mesh} & $s$ & Cells ($10^6$) & $r$ & Faces across bores & $R$ & '
           r'Error (\%) & GCI (\%) & $R_\infty$ \\', r'\midrule']
    notes = []
    design = None
    for r in rows:
        first = r['design'] != design
        if first and design is not None:
            out.append(r'\addlinespace')
        design = r['design']
        name = tex(design) if first else ''
        if r['note'].strip():
            notes.append(r'\item[' + tex(design) + '] ' + tex_note(r['note']) + '.')
        out.append(' & '.join([name, r['mesh'], r['s'], r['cells_M'], r.get('r_from_coarser', ''),
                               tex(r['faces_across_bores']).replace('; ', ', '), r['R_U0.3'],
                               r['medium_error_pct'], r['GCI_fine_pct'], r['R_best']]) + r' \\')
    out += [r'\bottomrule', r'\end{tabular}']
    if notes:
        out += [r'\begin{description}\footnotesize\setlength{\itemsep}{0pt}'] + notes + \
               [r'\end{description}']
    out += [r'\end{table}']
    return '\n'.join(out) + '\n'


def table_s4():
    rows = read('S4_runs.csv')
    # Misalignment runs are not part of this paper's run list.
    rows = [r for r in rows if r['purpose'].strip() != 'misalignment']
    out = [r'\begin{longtable}{@{}' + RAG + r'p{0.11\textwidth}' + RAG +
           r'p{0.21\textwidth}' + RAG + r'p{0.11\textwidth}rrrrrrr@{}}',
           r'\caption{Every run whose result the paper reports. Unless stated, coupled URANS with'
           r' SST $k$--$\omega$ $\gamma$ on the half model: 1500 steady iterations, a ramp of 100 steps'
           r' at $\Delta t/20$ and 10 steps each at $\Delta t/10$, $\Delta t/5$, $2\Delta t/5$ and $4\Delta t/5$'
           r' with 20 inner iterations, then 10 inner iterations per step, inflow turbulence intensity'
           r' \SI{5}{\percent} and viscosity ratio 10; $R$ averaged over the second half of the run.'
           r' s.d.\ and drift: standard deviation of $R$ over the averaging window, and the mean of its'
           r' second half less that of its first, in per cent of $R$. Imbalance: inlet-mouth less'
           r' outlet-mouth flow, relative to the mean inflow; $^\dagger$ above the $10^{-4}$ criterion'
           r' (scatter between inner-iteration reports, not a mean imbalance). Case suffixes: sst, plain SST;'
           r' rke, realizable $k$--$\varepsilon$; sst-gamma, the standard SST $k$--$\omega$ $\gamma$.'
           r'}\label{tab:S5}\\',
           r'\toprule', r'\textbf{Purpose} & \textbf{Case} & \textbf{Domain} & $U$ & $\Delta t$ (s) & '
           r'$t$ (s) & $R$ & s.d.\ (\%) & Drift (\%) & Imbalance \\', r'\midrule', r'\endfirsthead',
           r'\toprule', r'\textbf{Purpose} & \textbf{Case} & \textbf{Domain} & $U$ & $\Delta t$ (s) & '
           r'$t$ (s) & $R$ & s.d.\ (\%) & Drift (\%) & Imbalance \\', r'\midrule', r'\endhead',
           r'\bottomrule', r'\endlastfoot']
    purpose = None
    for r in rows:
        p = tex(r['purpose']) if r['purpose'] != purpose else ''
        if r['purpose'] != purpose and purpose is not None:
            out.append(r'\addlinespace')
        purpose = r['purpose']
        case = r'\texttt{' + tex(r['case']).replace(r'\_', r'\_\allowbreak{}') + '}'
        out.append(' & '.join([p, case, tex(r['domain']), r['U'], r['dt_s'], r['t_end_s'],
                               r['R'], r['sd_pct'], r['drift_pct'],
                               r'\num{' + r['imbalance'] + '}' +
                               (r'$^\dagger$' if r.get('imbalance_flag', '').strip() else '')]) + r' \\')
    out.append(r'\end{longtable}')
    return '\n'.join(out) + '\n'


if __name__ == '__main__':
    for name, fn in [('tab_S1.tex', table_s1), ('tab_S2.tex', table_s2), ('tab_S3.tex', table_s3),
                     ('tab_S4.tex', table_s4)]:
        with open(os.path.join(HERE, name), 'w', encoding='utf-8', newline='\n') as f:
            f.write('% Generated by make_supp_tables.py from cfd/results/supplemental/; do not edit.\n')
            f.write(fn())
        print('wrote', name)
