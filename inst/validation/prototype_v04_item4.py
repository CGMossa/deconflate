"""Prototype for v0.4 item 4, step 2 (and the new association rows of step 1).
Run from this directory: python3 prototype_v04_item4.py   (numpy, scipy)

Proposed input layout (one row = one estimate):
  diseases.csv     : disease, measure (prevalence | probability | incidence_rate), value
  associations.csv : disease, associated_with, among, adjusted_for, measure, value
                     disease / associated_with: a disease or a combination ("MAS;SCK")
                     among: a subgroup ("SCK", "no SCK", "SCK;no MET")
                     measure: OR, RR, RD, phi, prevalence, OR_ratio
  impacts.csv      : disease, among, adjusted_for, value (+ measure for event impacts)

A. Associations. Every row is turned into the probability of one event of
   disease combinations (e.g. P(LAM and MAS and SCK)), given the current fit
   where the conversion needs other probabilities (a fixed point). The joint
   distribution is the maximum-entropy one with all these event
   probabilities ("raking on events"), which for pairs only is the current
   IPF fit.
   Checks:
   A1. pairs only: raking on events = the current IPF fit (five diseases; 2024).
   A2. each kind of row recovers a known joint distribution exactly.
   A3. the current three-way term (a log-linear parameter, conditional on all
       other diseases) is not the stratum odds ratio ratio of a study.
   A4. conflicting rows are detected.
B. Additive impacts with combinations and subgroups. Unknowns: each
   disease's own impact and an interaction for each combination listed as a
   `disease`. Each row is one linear equation (crude contrast or linear
   projection, within `among`).
   Checks: B1 recovery of known impacts and interactions; B2 identification
   check; B3 "only vs healthy" vs "adjusted for all"; B4 Shapley split of
   n-way interactions; B5 a combination row on the 2024 inputs.
C. Event impacts with combinations (interaction terms on the log-hazard
   scale): recovery from crude hazard ratios of single diseases and a
   combination.
"""
import itertools
import time
import numpy as np
from scipy.optimize import brentq, fsolve
from python_reference import or_to_joint

np.set_printoptions(precision=6, linewidth=150, suppress=True)


# ---------------------------------------------------------------- helpers
def cells_of(n):
    return np.array(list(itertools.product([0, 1], repeat=n)), float)


class Pop:
    def __init__(self, ids):
        self.ids = list(ids)
        self.ix = {d: i for i, d in enumerate(ids)}
        self.cells = cells_of(len(ids))

    def ev(self, pos=(), neg=()):
        """Indicator of the event: all of `pos` present, all of `neg` absent."""
        m = np.ones(len(self.cells), bool)
        for d in pos:
            m &= self.cells[:, self.ix[d]] == 1
        for d in neg:
            m &= self.cells[:, self.ix[d]] == 0
        return m


def split(s):
    return [x.strip() for x in s.split(";")] if s else []


def among_sets(among):
    pos, neg = [], []
    for a in split(among):
        if a.startswith("no "):
            neg.append(a[3:].strip())
        else:
            pos.append(a)
    return pos, neg


def to_joint(meas, v, p1, p2):
    """P(1 and 2) from a 2x2 measure; p1 = P(disease), p2 = P(associated_with)."""
    if meas == "OR":
        return or_to_joint(p1, p2, v)
    if meas == "RR":       # P(1 | 2) / P(1 | not 2)
        return v * p1 / (1 - p2 + v * p2) * p2
    if meas == "RD":       # P(1 | 2) - P(1 | not 2)
        return (v * (1 - p2) + p1) * p2
    if meas == "phi":
        return p1 * p2 + v * np.sqrt(p1 * (1 - p1) * p2 * (1 - p2))
    raise ValueError(meas)


def odds_ratio(p11, p1, p2):
    return p11 * (1 - p1 - p2 + p11) / ((p1 - p11) * (p2 - p11))


def stratum_or(pop, p, g, h, apos, aneg):
    A = pop.ev(apos, aneg)
    pa = p[A].sum()
    pg = p[A & pop.ev(g)].sum() / pa
    ph = p[A & pop.ev(h)].sum() / pa
    pgh = p[A & pop.ev(g + h)].sum() / pa
    return odds_ratio(pgh, pg, ph)


def logistic_coef(pop, p, y_set, x_sets, apos=(), aneg=()):
    """Population logistic projection of 1{y_set} on [1, 1{x} for x in x_sets]
    within the subgroup; returns the coefficients (IRLS)."""
    A = pop.ev(apos, aneg)
    w = p[A]
    Z = np.column_stack([np.ones(A.sum())] + [pop.ev(x)[A].astype(float) for x in x_sets])
    yv = pop.ev(y_set)[A].astype(float)
    beta = np.zeros(Z.shape[1])
    for _ in range(200):
        mu = 1 / (1 + np.exp(-Z @ beta))
        g = Z.T @ (w * (yv - mu))
        H = (Z * (w * mu * (1 - mu))[:, None]).T @ Z
        step = np.linalg.solve(H, g)
        beta += step
        if np.abs(step).max() < 1e-13:
            break
    return beta


# ---------------------------------------------------------------- A. associations
def row_event(pop, p, r, prev):
    """The event (pos, neg) and its target probability for one row, given the
    current joint p (used where the conversion needs other probabilities)."""
    g, h = split(r["disease"]), split(r.get("with", ""))
    apos, aneg = among_sets(r.get("among", ""))
    m, v = r.get("measure", "OR"), r["value"]
    A = pop.ev(apos, aneg)
    pA = p[A].sum()
    if m == "prevalence":                        # P(g | h, among)
        HA = A & pop.ev(h)
        return (g + h + apos, aneg), v * p[HA].sum()
    pg = p[A & pop.ev(g)].sum() / pA             # margins within the subgroup
    ph = p[A & pop.ev(h)].sum() / pA
    if r.get("adjusted_for"):
        # Mapped to the stratum OR it implies at the current fit (fixed point):
        # target log OR = log v + (current log OR - current adjusted log OR).
        adj = [[a] for a in pop.ids if a not in g + h] if r["adjusted_for"] == "all" else \
            [[a] for a in split(r["adjusted_for"])]
        cur_adj = logistic_coef(pop, p, g, [h] + adj, apos, aneg)[1]
        cur = np.log(stratum_or(pop, p, g, h, apos, aneg))
        v = np.exp(np.log(v) + cur - cur_adj)
        m = "OR"
    if m == "OR_ratio":                          # OR among `among` / OR among its complement
        assert len(apos) == 1 and not aneg, "OR_ratio needs one disease in `among`"
        # Both stratum ORs are functions of t = P(g, h, S) once the margins
        # P(S), P(g, S), P(h, S), P(g), P(h) and P(g, h) are fixed (current fit):
        # solve ratio(t) = v for t.
        S = apos
        pS = p[pop.ev(S)].sum()
        pgS, phS = p[pop.ev(g + S)].sum(), p[pop.ev(h + S)].sum()
        pgT, phT, pghT = p[pop.ev(g)].sum(), p[pop.ev(h)].sum(), p[pop.ev(g + h)].sum()
        def ratio(t):
            o1 = odds_ratio(t / pS, pgS / pS, phS / pS)
            o0 = odds_ratio((pghT - t) / (1 - pS), (pgT - pgS) / (1 - pS), (phT - phS) / (1 - pS))
            return np.log(o1) - np.log(o0) - np.log(v)
        lo = max(0.0, pgS + phS - pS, pghT - min(pgT - pgS, phT - phS)) + 1e-12
        hi = min(pgS, phS, pghT, pghT - (pgT - pgS + phT - phS - (1 - pS))) - 1e-12
        t = brentq(ratio, lo, hi, xtol=1e-15)
        return (g + h + S, []), t
    q = to_joint(m, v, pg, ph)                   # P(g and h | among)
    return (g + h + apos, aneg), q * pA


def fit(pop, prevalences, rows, start=None, tol=1e-12, max_outer=500):
    """Maximum-entropy joint with the given prevalences and association rows."""
    P = np.array([prevalences[d] for d in pop.ids])
    if start is None:
        p = np.prod(np.where(pop.cells == 1, P, 1 - P), axis=1)
    else:
        p = start.copy()
    single = [((d,), ()) for d in pop.ids]
    single_t = list(P)
    masks_single = [pop.ev(e[0], e[1]) for e in single]
    for outer in range(max_outer):
        events = [row_event(pop, p, r, prevalences) for r in rows]
        masks = masks_single + [pop.ev(e[0], e[1]) for e, _ in events]
        targets = single_t + [t for _, t in events]
        for sweep in range(200):
            for m, t in zip(masks, targets):
                s = p[m].sum()
                p[m] *= t / s
                p[~m] *= (1 - t) / (1 - s)
            err = max(abs(p[m].sum() - t) for m, t in zip(masks, targets))
            if err < tol:
                break
        # outer convergence: the conversions no longer move
        new = [row_event(pop, p, r, prevalences)[1] for r in rows]
        if max([abs(a - b) for a, b in zip(new, [t for _, t in events])] + [0]) < tol and err < tol:
            return p, outer + 1, err
    return p, max_outer, err


def row_value(pop, p, r):
    """The value a row describes, evaluated on a joint distribution p."""
    g, h = split(r["disease"]), split(r.get("with", ""))
    apos, aneg = among_sets(r.get("among", ""))
    m = r.get("measure", "OR")
    A = pop.ev(apos, aneg)
    if m == "prevalence":
        return p[A & pop.ev(g + h)].sum() / p[A & pop.ev(h)].sum()
    if r.get("adjusted_for"):
        adj = [[a] for a in pop.ids if a not in g + h] if r["adjusted_for"] == "all" else \
            [[a] for a in split(r["adjusted_for"])]
        return np.exp(logistic_coef(pop, p, g, [h] + adj, apos, aneg)[1])
    if m == "OR_ratio":
        return stratum_or(pop, p, g, h, apos, []) / stratum_or(pop, p, g, h, [], apos)
    pA = p[A].sum()
    pg = p[A & pop.ev(g)].sum() / pA
    ph = p[A & pop.ev(h)].sum() / pA
    pgh = p[A & pop.ev(g + h)].sum() / pA
    if m == "OR":
        return odds_ratio(pgh, pg, ph)
    if m == "RR":
        return (pgh / ph) / ((pg - pgh) / (1 - ph))
    if m == "RD":
        return pgh / ph - (pg - pgh) / (1 - ph)
    if m == "phi":
        return (pgh - pg * ph) / np.sqrt(pg * (1 - pg) * ph * (1 - ph))


def loglinear_joint(pop, main, terms):
    """p(d) proportional to exp(sum main_i d_i + sum_S lambda_S prod_{i in S} d_i)."""
    lp = pop.cells @ np.array([main[d] for d in pop.ids])
    for S, lam in terms.items():
        lp = lp + lam * pop.ev(S.split(";"))
    p = np.exp(lp - lp.max())
    return p / p.sum()


print("=" * 78)
print("A. Association rows")
print("=" * 78)

# Five-disease population (as inst/extdata/five_diseases: every pair given).
ids5 = ["LAM", "MAS", "MET", "SCK", "RP"]
pop5 = Pop(ids5)
P5 = {"LAM": 0.25, "MAS": 1 - np.exp(-0.30), "MET": 0.10, "SCK": 0.35, "RP": 0.06}
rows5 = [dict(disease="LAM", **{"with": "MAS"}, measure="OR", value=1.8),
         dict(disease="LAM", **{"with": "SCK"}, measure="OR", value=2.0),
         dict(disease="MET", **{"with": "RP"}, measure="OR", value=6.0),
         dict(disease="MET", **{"with": "SCK"}, measure="RR", value=1.5),
         dict(disease="MAS", **{"with": "SCK"}, measure="prevalence", value=0.35),   # was cond_prob
         dict(disease="MAS", **{"with": "MET"}, measure="RD", value=0.08),
         dict(disease="LAM", **{"with": "MET"}, measure="phi", value=0.05),
         dict(disease="RP", **{"with": "SCK"}, measure="OR", value=30 * 610 / (40 * 320)),
         dict(disease="LAM", **{"with": "RP"}, measure="OR", value=1.0),
         dict(disease="MAS", **{"with": "RP"}, measure="OR", value=1.5)]

# A1. Pairs only: equals the current IPF (python reference of the package).
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):
    import reference_five_diseases as f5
p_old = f5.ipf(f5.P, f5.J, three=(("LAM", "MAS", "SCK"), 1.0), iters=4000)
t0 = time.time()
p_new, n_outer, err = fit(pop5, P5, rows5)
print(f"A1. five diseases, pairs only: max |new - current IPF| = {np.abs(p_new - p_old).max():.2e}"
      f" (outer {n_outer}, residual {err:.1e}, {time.time() - t0:.2f} s)")
# Note: the `cond_prob` row P(MAS | SCK) = 0.35 is the `prevalence` row MAS with SCK.

# 2024 global dairy (12 diseases, 38 pairs, 28 unknown): yield aggregate.
import reference_2024_analysis as r24   # (prints its own values; ignored)
D24, P24 = r24.D, r24.P
pop24 = Pop(D24)
rows24 = [dict(disease=D24[i], **{"with": D24[k]}, measure="OR", value=r24.central(spec))
          for i, k, spec in r24.pairs]
t0 = time.time()
p24, n_outer, err = fit(pop24, dict(zip(D24, P24)), rows24, tol=1e-11)
print(f"\nA1. 2024, 12 diseases, 38 pairs + 28 unknown: outer {n_outer}, residual {err:.1e}, {time.time() - t0:.1f} s")


def crude_A(pop, p):
    n = len(pop.ids)
    Pm = pop.cells.T @ p
    J = (pop.cells * p[:, None]).T @ pop.cells
    S = J - np.outer(Pm, Pm)
    np.fill_diagonal(S, Pm * (1 - Pm))
    return S / np.diag(S)[:, None], Pm


A24, Pm24 = crude_A(pop24, p24)
y24 = np.array([r24.central(r24.YIELD[d]) for d in D24])
b24 = np.linalg.solve(A24, y24)
print(f"    yield aggregate {Pm24 @ b24:.10f} (reference_v040.py: 6.9150664186)")

# A2. Each kind of row recovers a known joint exactly. Truth: the pairwise
# log-linear model of the five diseases plus a LAM:MAS:SCK term.
lam = {"LAM;MAS": 0.6, "LAM;SCK": 0.7, "MET;RP": 1.8, "MET;SCK": 0.4, "MAS;SCK": 0.3,
       "MAS;MET": 0.5, "LAM;MET": 0.2, "RP;SCK": 0.35, "MAS;RP": 0.4, "LAM;MAS;SCK": 0.5}
main = {"LAM": -1.3, "MAS": -1.2, "MET": -2.3, "SCK": -0.8, "RP": -3.0}
p_true = loglinear_joint(pop5, main, lam)
P_true = {d: p_true[pop5.ev([d])].sum() for d in ids5}
pairs = [r.split(";") for r in lam if r.count(";") == 1]
pair_rows = [dict(disease=a, **{"with": b}, measure="OR",
                  value=row_value(pop5, p_true, dict(disease=a, **{"with": b}, measure="OR")))
             for a, b in pairs]
print("\nA2. Truth: pairwise log-linear model + a LAM:MAS:SCK term (lambda 0.5).")
print("    Prevalences and the 9 pair ORs from the truth, plus ONE row of each kind below,")
print("    each with its value computed from the truth:")
p_pairs_only, _, _ = fit(pop5, P_true, pair_rows)
print(f"    {'(pairs only: no three-way information)':52s} max |fit - truth| = {np.abs(p_pairs_only - p_true).max():.2e}")
kinds = [
    ("LAM by MAS;SCK (OR, combination)", dict(disease="LAM", **{"with": "MAS;SCK"}, measure="OR")),
    ("LAM by MAS among SCK (OR, stratum)", dict(disease="LAM", **{"with": "MAS"}, among="SCK", measure="OR")),
    ("LAM by MAS among no SCK (OR, stratum)", dict(disease="LAM", **{"with": "MAS"}, among="no SCK", measure="OR")),
    ("LAM by MAS among SCK (OR_ratio)", dict(disease="LAM", **{"with": "MAS"}, among="SCK", measure="OR_ratio")),
    ("LAM by MAS among SCK (RR, stratum)", dict(disease="LAM", **{"with": "MAS"}, among="SCK", measure="RR")),
    ("LAM by MAS;SCK (RD, combination)", dict(disease="LAM", **{"with": "MAS;SCK"}, measure="RD")),
    ("LAM with MAS;SCK (prevalence: P(LAM | MAS, SCK))", dict(disease="LAM", **{"with": "MAS;SCK"}, measure="prevalence")),
    ("LAM;MAS;SCK (prevalence of the combination)", dict(disease="LAM;MAS;SCK", measure="prevalence")),
]
for lab, r in kinds:
    r = dict(r, value=row_value(pop5, p_true, r))
    pf, no, er = fit(pop5, P_true, pair_rows + [r])
    print(f"    {lab:52s} value {r['value']:.6f}  max |fit - truth| = {np.abs(pf - p_true).max():.2e}  (outer {no})")

# A3. The current three-way term vs a study's stratum OR ratio.
p3 = f5.ipf(f5.P, f5.J, three=(("LAM", "MAS", "SCK"), 1.5), iters=4000)
r_ratio = dict(disease="LAM", **{"with": "MAS"}, among="SCK", measure="OR_ratio")
print("\nA3. The current three-way term (exp(lambda_LMS) = 1.5, conditional on all other diseases):")
print(f"    the stratum OR ratio a study would report: {row_value(pop5, p3, r_ratio):.6f}")
pr, _, _ = fit(pop5, P5, rows5 + [dict(r_ratio, value=1.5)])
cells_lms = pop5.ev(["LAM", "MAS", "SCK"])
print(f"    with OR_ratio = 1.5 as a row: P(LAM, MAS, SCK) = {pr[cells_lms].sum():.6f} "
      f"(current term: {p3[cells_lms].sum():.6f}; none: {p_new[cells_lms].sum():.6f})")

# A4. Conflicting rows: both stratum ORs and an OR_ratio that does not match.
r_s = dict(disease="LAM", **{"with": "MAS"}, among="SCK", measure="OR", value=2.4)
r_n = dict(disease="LAM", **{"with": "MAS"}, among="no SCK", measure="OR", value=1.6)
ok = fit(pop5, P5, [r for r in rows5 if not (r["disease"] == "LAM" and r["with"] == "MAS")] + [r_s, r_n])
print("\nA4. Stratum ORs 2.4 (SCK) and 1.6 (no SCK) instead of the crude LAM-MAS OR:")
print(f"    fitted: {row_value(pop5, ok[0], r_s):.6f}, {row_value(pop5, ok[0], r_n):.6f}; implied crude OR "
      f"{row_value(pop5, ok[0], dict(disease='LAM', **{'with': 'MAS'}, measure='OR')):.6f}; outer {ok[1]}, residual {ok[2]:.1e}")
bad = fit(pop5, P5, [r for r in rows5 if not (r["disease"] == "LAM" and r["with"] == "MAS")] +
          [r_s, r_n, dict(disease="LAM", **{"with": "MAS"}, among="SCK", measure="OR_ratio", value=2.0)], max_outer=300)
print(f"    plus OR_ratio 2.0 (2.4 / 1.6 = 1.5 conflicts): outer {bad[1]} (limit 300), residual {bad[2]:.1e},"
      f" rows reproduced: {[round(row_value(pop5, bad[0], r), 3) for r in (r_s, r_n)]}")

# ---------------------------------------------------------------- B. impacts
print("\n" + "=" * 78)
print("B. Additive impacts with combinations and subgroups")
print("=" * 78)


def impact_row_coefs(pop, p, r, unknowns):
    """Coefficients of one impact row on the unknowns (disease ids and
    combinations 'A;B'): the row's contrast applied to each unknown's
    indicator variable."""
    g = split(r["disease"])
    apos, aneg = among_sets(r.get("among", ""))
    A = pop.ev(apos, aneg)
    w = p[A] / p[A].sum()
    G = pop.ev(g)[A].astype(float)
    if r.get("adjusted_for"):
        adj = [x for x in pop.ids if x not in g] if r["adjusted_for"] == "all" else split(r["adjusted_for"])
        Z = np.column_stack([np.ones(A.sum()), G] + [pop.ev([a])[A].astype(float) for a in adj])
        ZtW = (Z * w[:, None]).T
        M = np.linalg.solve(ZtW @ Z, ZtW)                       # projection coefficients
        return np.array([(M @ pop.ev(split(u))[A].astype(float))[1] for u in unknowns])
    p1, p0 = w * G, w * (1 - G)
    return np.array([(p1 @ pop.ev(split(u))[A]) / p1.sum() - (p0 @ pop.ev(split(u))[A]) / p0.sum()
                     for u in unknowns])


def outcome_mean_contrast(pop, p, r, b, delta):
    """The value an impact row describes under the additive model with impacts
    b and interactions delta (computed directly on the cells)."""
    y = pop.cells @ np.array([b[d] for d in pop.ids])
    for S, v in delta.items():
        y = y + v * pop.ev(split(S))
    g = split(r["disease"])
    apos, aneg = among_sets(r.get("among", ""))
    A = pop.ev(apos, aneg)
    w = p[A] / p[A].sum()
    G = pop.ev(g)[A].astype(float)
    if r.get("adjusted_for"):
        adj = [x for x in pop.ids if x not in g] if r["adjusted_for"] == "all" else split(r["adjusted_for"])
        Z = np.column_stack([np.ones(A.sum()), G] + [pop.ev([a])[A].astype(float) for a in adj])
        ZtW = (Z * w[:, None]).T
        return np.linalg.solve(ZtW @ Z, ZtW @ y[A])[1]
    return (w * G) @ y[A] / (w * G).sum() - (w * (1 - G)) @ y[A] / (w * (1 - G)).sum()


def solve_impacts(pop, p, rows):
    unknowns = list(pop.ids) + [r["disease"] for r in rows if ";" in r["disease"]]
    unknowns = list(dict.fromkeys(unknowns))
    Am = np.array([impact_row_coefs(pop, p, r, unknowns) for r in rows])
    y = np.array([r["value"] for r in rows])
    rank = np.linalg.matrix_rank(Am, tol=1e-9)
    if rank < len(unknowns):
        # Which unknowns are not identified: those with a component in the null space.
        _, s, Vt = np.linalg.svd(Am)
        null = Vt[rank:]
        bad = [u for u, c in zip(unknowns, np.abs(null).max(0)) if c > 1e-8]
        return unknowns, None, rank, bad
    theta = np.linalg.lstsq(Am, y, rcond=None)[0]
    return unknowns, theta, rank, []


# Truth: impacts and interactions on the five-disease joint (pairs given, no
# three-way association), interactions LAM:MAS = 2 and LAM:MAS:SCK = -1.
p5 = p_new
b_true = {"LAM": 4.8, "MAS": 3.0, "MET": 5.0, "SCK": 1.5, "RP": 2.5}
d_true = {"LAM;MAS": 2.0, "LAM;MAS;SCK": -1.0}
rowsB = [dict(disease=d) for d in ids5] + [dict(disease="LAM;MAS"), dict(disease="LAM;MAS;SCK", adjusted_for="all")]
for r in rowsB:
    r["value"] = outcome_mean_contrast(pop5, p5, r, b_true, d_true)
unk, theta, rank, bad = solve_impacts(pop5, p5, rowsB)
print("B1. Truth: impacts", b_true, "interactions", d_true)
print("    Rows: crude for each disease, crude LAM;MAS, LAM;MAS;SCK adjusted for all; values from the truth:")
for r in rowsB:
    print(f"      {r['disease']:12s} {('adjusted for ' + r['adjusted_for']) if r.get('adjusted_for') else 'crude':18s} {r['value']:.6f}")
print("    recovered:", dict(zip(unk, np.round(theta, 10))), f"rank {rank}/{len(unk)}")
# Other ways to inform the interaction: a stratum row and an adjusted combination row.
rowsB2 = [dict(disease=d) for d in ids5] + [dict(disease="LAM;MAS", adjusted_for="all"),
                                             dict(disease="LAM", among="MAS"),
                                             dict(disease="LAM;MAS;SCK")]
for r in rowsB2:
    r["value"] = outcome_mean_contrast(pop5, p5, r, b_true, d_true)
unk2, theta2, rank2, _ = solve_impacts(pop5, p5, rowsB2)
print("    Rows: crude singles, LAM;MAS adjusted for all, LAM among MAS (crude), LAM;MAS;SCK crude (8 rows, 7 unknowns):")
print("    recovered:", dict(zip(unk2, np.round(theta2, 10))), f"rank {rank2}/{len(unk2)}")

# B2. Identification: a disease without a row of its own.
rowsB3 = [dict(disease=d, value=1.0) for d in ids5 if d != "LAM"] + [dict(disease="LAM;MAS", value=1.0)]
unk3, theta3, rank3, bad3 = solve_impacts(pop5, p5, rowsB3)
print(f"\nB2. Rows for MAS, MET, SCK, RP and LAM;MAS only: rank {rank3}/{len(unk3)}; not identified: {bad3}")
rowsB4 = [dict(disease=d, value=1.0) for d in ids5] + [dict(disease="LAM;MAS", value=1.0),
                                                        dict(disease="LAM;MAS", among="SCK", value=1.0)]
unk4, theta4, rank4, bad4 = solve_impacts(pop5, p5, rowsB4)
print(f"    Crude singles + LAM;MAS crude + LAM;MAS among SCK: rank {rank4}/{len(unk4)} (7 rows, 6 unknowns: over-determined, item 4 step 3)")

# B3. "CM only vs healthy" vs "adjusted for all" with an interaction present.
y = pop5.cells @ np.array([b_true[d] for d in ids5]) + 2.0 * pop5.ev(["LAM", "MAS"])
only = pop5.ev(["LAM"], [d for d in ids5 if d != "LAM"])
healthy = pop5.ev([], ids5)
print(f"\nB3. LAM only vs healthy: {y[only].mean() - y[healthy].mean():.6f} (= LAM's own impact {b_true['LAM']}); "
      f"LAM adjusted for all: {outcome_mean_contrast(pop5, p5, dict(disease='LAM', adjusted_for='all'), b_true, {'LAM;MAS': 2.0}):.6f}"
      f" (interaction LAM:MAS = 2)")

# B4. Shapley split of n-way interactions: closed form vs permutations.
Pm5 = pop5.cells.T @ p5
closed = {d: Pm5[pop5.ix[d]] * b_true[d] for d in ids5}
for S, v in d_true.items():
    pS = p5[pop5.ev(split(S))].sum()
    for d in split(S):
        closed[d] += v * pS / len(split(S))
import math
perm = {d: 0.0 for d in ids5}
for c, pc in zip(pop5.cells, p5):
    present = [d for d in ids5 if c[pop5.ix[d]] == 1]
    def loss(subset):
        v = sum(b_true[d] for d in subset)
        for S, dv in d_true.items():
            if set(split(S)) <= set(subset):
                v += dv
        return v
    nf = math.factorial(len(present))
    for order in itertools.permutations(present):
        acc = []
        for d in order:
            before = loss(acc)
            acc.append(d)
            perm[d] += pc * (loss(acc) - before) / nf
print("\nB4. Shapley contributions (closed form: each interaction split equally among its diseases):")
print("    closed form ", {d: round(v, 10) for d, v in closed.items()})
print("    permutations", {d: round(v, 10) for d, v in perm.items()})
print(f"    aggregate {sum(closed.values()):.10f} = E[Y] {p5 @ y + (-1.0) * p5 @ pop5.ev(['LAM', 'MAS', 'SCK']):.10f}")

# B5. 2024 yield with a (hypothetical) combination row SCK;SCM.
print("\nB5. 2024 yield (unknown pairs, global fit) with a hypothetical crude row for cows with both SCK and SCM:")
rows24i = [dict(disease=d, value=v) for d, v in zip(D24, y24)]
base = solve_impacts(pop24, p24, rows24i)
pm = pop24.cells.T @ p24
pSS = p24[pop24.ev(["SCK", "SCM"])].sum()
print(f"    without it: aggregate {pm @ base[1]:.6f}; P(SCK and SCM) = {pSS:.6f}")
for v in (12.0, 14.69, 18.0):
    u, th, rk, _ = solve_impacts(pop24, p24, rows24i + [dict(disease="SCK;SCM", value=v)])
    agg = pm @ th[:len(D24)] + th[-1] * pSS
    print(f"    SCK;SCM crude = {v:5.2f}: interaction {th[-1]:+.4f}, SCK {th[D24.index('SCK')]:.4f}, "
          f"SCM {th[D24.index('SCM')]:.4f}, aggregate {agg:.6f}")
# The value of the row with no interaction: what the model implies.
r_imp = dict(disease="SCK;SCM")
implied = impact_row_coefs(pop24, p24, r_imp, list(D24)) @ base[1]
print(f"    (the value the inputs imply with no interaction: {implied:.6f})")

# ---------------------------------------------------------------- C. event impacts
print("\n" + "=" * 78)
print("C. Event impacts with a combination (interaction on the log-hazard scale)")
print("=" * 78)


def event_hr(pop, p, r, beta, gamma):
    """Crude hazard ratio of a row: ratio of mean hazard multipliers among
    animals with and without the row's disease or combination (snapshot)."""
    lm = pop.cells @ np.array([beta[d] for d in pop.ids])
    for S, v in gamma.items():
        lm = lm + v * pop.ev(split(S))
    mult = np.exp(lm)
    G = pop.ev(split(r["disease"]))
    return (p[G] @ mult[G] / p[G].sum()) / (p[~G] @ mult[~G] / p[~G].sum())


beta_true = {"LAM": np.log(1.6), "MAS": np.log(1.5), "MET": np.log(1.1), "SCK": np.log(1.2), "RP": np.log(1.4)}
gamma_true = {"LAM;MAS": np.log(1.3)}
rowsC = [dict(disease=d) for d in ids5] + [dict(disease="LAM;MAS")]
for r in rowsC:
    r["value"] = event_hr(pop5, p5, r, beta_true, gamma_true)
names = list(ids5) + ["LAM;MAS"]


def F(x):
    beta = dict(zip(ids5, x[:5]))
    gamma = {"LAM;MAS": x[5]}
    return [np.log(event_hr(pop5, p5, r, beta, gamma)) - np.log(r["value"]) for r in rowsC]


x = fsolve(F, np.zeros(6), xtol=1e-13)
print("Truth: hazard multipliers", {d: round(np.exp(v), 4) for d, v in beta_true.items()}, "and LAM:MAS x1.3")
print("Rows: crude HRs", {r["disease"]: round(r["value"], 6) for r in rowsC})
print("Recovered:", dict(zip(names, np.round(np.exp(x), 8))), f"max residual {np.abs(F(x)).max():.1e}")
