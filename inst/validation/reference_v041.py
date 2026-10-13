"""Reference values for deconflate 0.4, item 4 (one model, the new input
layout, rows about combinations and subgroups, several rows per item).

Independent implementation in Python (numpy, scipy) of the definitions in
?cm_associations, ?cm_impacts and ?deconflate. Run from this directory:

    python3 reference_v041.py

It reads the five-disease example (inst/extdata/five_diseases) and prints
the values checked by tests/testthat/test-import.R and test-event.R (the
self-consistency checks of the definitions are in test-item4.R).

Definitions (all vs none: a row about a combination compares animals with
all of its diseases to animals with none of them):
  associations  G = disease, H = associated_with, A = among
    OR/RR/RD/phi  2x2 table of all G by (all H vs none H) within A
    prevalence    P(all G | all H, A)
    OR_ratio      OR(G, H | S) / OR(G, H | not S)
    adjusted OR   exp(coef of H) in the logistic regression of G on H and the
                  adjustment diseases, within A
    Rows fixing the same event are pooled (DerSimonian-Laird on the log OR of
    their table; a row without a distribution is exact).
  additive impacts  y(c) = sum b_d c_d + sum delta_T all_T(c)
    crude    E[y | all G, A] - E[y | none G, A]
    adjusted coefficient of all_G in the linear regression within A on all_G,
             (combination) its members and smaller combinations, and the
             adjustment diseases; "all" = every other disease and term.
    several rows: weighted least squares (1/SD^2), exact rows as constraints
  event impacts  log hazard multiplier sum beta_d c_d + sum gamma_T all_T(c)
    crude    all G vs none G within A (HR: mean multipliers; RR/OR/RD: risks)
    adjusted within strata of the adjustment diseases, the interaction
             contrast over the subsets of G, weights 1 / sum_U 1/d_U;
             HR adjusted for all = the term's own value
"""
import csv
import itertools
import os
import numpy as np
from scipy import stats
from scipy.optimize import brentq, least_squares

np.set_printoptions(precision=10, linewidth=150)
HERE = os.path.dirname(os.path.abspath(__file__))
FIVE = os.path.join(HERE, "..", "extdata", "five_diseases")


# ------------------------------------------------------------------ inputs
def read(name):
    with open(os.path.join(FIVE, name), newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    for i, r in enumerate(rows):
        r["row"] = i + 1
    return rows


def num(x):
    return None if x is None or x.strip() == "" else float(x)


def quantile_fn(r):
    """Quantile function of a row's distribution (None: exact)."""
    d = (r.get("dist") or "").strip()
    p = [num(r.get(k, "")) for k in ("p1", "p2", "p3", "p4")]
    if d in ("", "fixed"):
        return None
    if d == "normal":
        lo = -np.inf if p[2] is None else p[2]
        hi = np.inf if p[3] is None else p[3]
        a, b = (lo - p[0]) / p[1], (hi - p[0]) / p[1]
        return lambda u: stats.truncnorm.ppf(u, a, b, loc=p[0], scale=p[1])
    if d == "lognormal":
        return lambda u: np.exp(stats.norm.ppf(u, p[0], p[1]))
    if d == "lognormal_ci":
        lev = 0.95 if p[3] is None else p[3]
        z = stats.norm.ppf(1 - (1 - lev) / 2)
        s = (np.log(p[2]) - np.log(p[1])) / (2 * z)
        return lambda u: np.exp(stats.norm.ppf(u, np.log(p[0]), s))
    if d in ("pert", "pert_mean"):
        lam = 4 if p[3] is None else p[3]
        mn, mx = p[0], p[2]
        mode = p[1] if d == "pert" else ((lam + 2) * p[1] - mn - mx) / lam
        a = 1 + lam * (mode - mn) / (mx - mn)
        b = 1 + lam * (mx - mode) / (mx - mn)
        return lambda u: mn + (mx - mn) * stats.beta.ppf(u, a, b)
    if d == "beta":
        mn = 0 if p[2] is None else p[2]
        mx = 1 if p[3] is None else p[3]
        return lambda u: mn + (mx - mn) * stats.beta.ppf(u, p[0], p[1])
    if d == "uniform":
        return lambda u: p[0] + (p[1] - p[0]) * u
    raise ValueError(d)


def se_on(qf, z):
    """SD on a scale z: half the distance between the transformed 15.9% and
    84.1% quantiles."""
    if qf is None:
        return 0.0
    return abs(z(qf(stats.norm.cdf(1))) - z(qf(stats.norm.cdf(-1)))) / 2


def split(x):
    return [s.strip() for s in x.split(";") if s.strip()] if x else []


class Pop:
    def __init__(self, ids, P):
        self.ids = list(ids)
        self.P = np.array(P, float)
        self.cells = np.array(list(itertools.product([0, 1], repeat=len(ids))), float)

    def col(self, d):
        return self.cells[:, self.ids.index(d)] == 1

    def allm(self, s):
        m = np.ones(len(self.cells), bool)
        for d in s:
            m &= self.col(d)
        return m

    def nonem(self, s):
        m = np.ones(len(self.cells), bool)
        for d in s:
            m &= ~self.col(d)
        return m

    def among(self, text, named):
        pos, neg = [], []
        others = False
        for it in split(text):
            if it.lower() == "no others":
                others = True
            elif it.lower().startswith("no "):
                neg.append(it[3:].strip())
            else:
                pos.append(it)
        if others:
            neg += [d for d in self.ids if d not in named + pos + neg]
        return pos, neg


def odds_ratio(p11, p1, p2):
    return p11 * (1 - p1 - p2 + p11) / ((p1 - p11) * (p2 - p11))


def or_to_joint(o, p1, p2):
    if abs(o - 1) < 1e-12:
        return p1 * p2
    lo, hi = max(0, p1 + p2 - 1), min(p1, p2)
    return brentq(lambda t: np.log(odds_ratio(t, p1, p2)) - np.log(o), lo + 1e-15, hi - 1e-15, xtol=1e-16)


def measure_joint(m, v, p1, p2):
    """P(G and H | R) from a 2x2 measure (G outcome, H comparison)."""
    if m == "OR":
        return or_to_joint(v, p1, p2)
    if m == "RR":
        return v * p1 / (1 - p2 + v * p2) * p2
    if m == "RD":
        return (v * (1 - p2) + p1) * p2
    if m == "phi":
        return p1 * p2 + v * np.sqrt(p1 * (1 - p1) * p2 * (1 - p2))
    if m == "prevalence":
        return v * p2
    raise ValueError(m)


# ------------------------------------------------------------ associations
def spec(pop, r):
    g, h = split(r["disease"]), split(r.get("associated_with", ""))
    apos, aneg = pop.among(r.get("among", ""), g + h)
    adj = split(r.get("adjusted_for", ""))
    if adj == ["all"]:
        adj = [d for d in pop.ids if d not in g + h + apos + aneg]
    else:
        adj = [d for d in adj if d in pop.ids]
    m = r.get("measure") or "OR"
    kind = "ratio" if m == "OR_ratio" else ("prevalence" if m == "prevalence" and not h else "table")
    pos = sorted(set(g + h + apos))
    return dict(row=r["row"], g=g, h=h, apos=apos, aneg=aneg, adj=adj, m=m, v=float(r["value"]),
                kind=kind, qf=quantile_fn(r), event=(tuple(pos), tuple(sorted(aneg))),
                table=(tuple(sorted(g + h)) if len(g) == 1 and len(h) == 1 else (tuple(g), tuple(h)),
                       tuple(apos), tuple(aneg)))


def table_probs(pop, s, p):
    A = pop.allm(s["apos"]) & pop.nonem(s["aneg"])
    aH = pop.allm(s["h"])
    R = A & (aH | pop.nonem(s["h"])) if len(s["h"]) > 1 else A
    pR = p[R].sum()
    return dict(pR=pR, pG=p[R & pop.allm(s["g"])].sum() / pR, pH=p[R & aH].sum() / pR,
                pA=p[A].sum(), E=pop.allm(s["g"]) & aH & A)


def logistic_coef(pop, p, s):
    A = pop.allm(s["apos"]) & pop.nonem(s["aneg"])
    w = p[A]
    Z = np.column_stack([np.ones(A.sum())] + [pop.col(x)[A].astype(float) for x in s["h"] + s["adj"]])
    y = pop.allm(s["g"])[A].astype(float)
    beta = np.zeros(Z.shape[1])
    for _ in range(100):
        mu = 1 / (1 + np.exp(-Z @ beta))
        step = np.linalg.solve((Z * (w * mu * (1 - mu))[:, None]).T @ Z, Z.T @ (w * (y - mu)))
        beta += step
        if np.abs(step).max() < 1e-13:
            break
    return beta[1]


def ratio_t(pop, s, v, p):
    S = s["apos"] + s["aneg"]
    neg = bool(s["aneg"])
    if neg:
        v = 1 / v
    g, h = s["g"], s["h"]
    pS = p[pop.allm(S)].sum()
    pgS, phS = p[pop.allm(g + S)].sum(), p[pop.allm(h + S)].sum()
    pgT, phT, pghT = p[pop.allm(g)].sum(), p[pop.allm(h)].sum(), p[pop.allm(g + h)].sum()
    a, b, mm = pgT - pgS, phT - phS, 1 - pS
    f = lambda t: (np.log(odds_ratio(t / pS, pgS / pS, phS / pS))
                   - np.log(odds_ratio((pghT - t) / mm, a / mm, b / mm)) - np.log(v))
    lo = max(0, pgS + phS - pS, pghT - a, pghT - b)
    hi = min(pgS, phS, pghT, pghT - (a + b - mm))
    eps = 1e-12 * (hi - lo)
    t = brentq(f, lo + eps, hi - eps, xtol=1e-16)
    return pghT - t if neg else t


def row_z(pop, s, v, p, scale):
    tp = table_probs(pop, s, p)
    if scale == "logOR":
        if s["adj"]:
            crude = np.log(odds_ratio(p[tp["E"]].sum() / tp["pR"], tp["pG"], tp["pH"]))
            return np.log(v) + crude - logistic_coef(pop, p, s)
        if s["m"] == "OR":
            return np.log(v)
        q = measure_joint(s["m"], v, tp["pG"], tp["pH"])
        return np.log(odds_ratio(q, tp["pG"], tp["pH"]))
    if scale == "logit_prev":
        return np.log(v / (1 - v))
    if scale == "log_ratio":
        return np.log(v)
    raise ValueError(scale)


def z_target(pop, s, z, p, scale):
    tp = table_probs(pop, s, p)
    if scale == "logOR":
        return or_to_joint(np.exp(z), tp["pG"], tp["pH"]) * tp["pR"]
    if scale == "logit_prev":
        return 1 / (1 + np.exp(-z)) * tp["pA"]
    if scale == "log_ratio":
        return ratio_t(pop, s, np.exp(z), p)


def pool(z, se, method="random"):
    z, se = np.array(z), np.array(se)
    if (se == 0).any():
        return z[np.argmax(se == 0)], dict(exact=True)
    if len(z) == 1:
        return z[0], dict(exact=False)
    w = 1 / se ** 2
    zf = (w * z).sum() / w.sum()
    Q = (w * (z - zf) ** 2).sum()
    df = len(z) - 1
    tau2 = max(0, (Q - df) / (w.sum() - (w ** 2).sum() / w.sum())) if method == "random" else 0
    ws = 1 / (se ** 2 + tau2)
    return (ws * z).sum() / ws.sum(), dict(Q=Q, df=df, tau2=tau2, I2=max(0, (Q - df) / Q) if Q > 0 else 0,
                                           se=np.sqrt(1 / ws.sum()))


def groups_of(specs):
    keys = []
    for s in specs:
        if s["event"] not in keys:
            keys.append(s["event"])
    out = []
    for k in keys:
        ss = [s for s in specs if s["event"] == k]
        kinds = {s["kind"] for s in ss}
        if kinds == {"table"} and len({s["table"] for s in ss}) == 1:
            scale = "logOR"
        elif kinds == {"prevalence"}:
            scale = "logit_prev"
        else:
            scale = "log_ratio"
        out.append(dict(event=k, specs=ss, scale=scale))
    return out


def fit(pop, specs, method="random", tol=1e-13, max_sweeps=20000):
    groups = groups_of(specs)
    p = np.prod(np.where(pop.cells == 1, pop.P, 1 - pop.P), axis=1)
    masks = [pop.allm(gr["event"][0]) & pop.nonem(gr["event"][1]) for gr in groups]
    margins = [pop.cells[:, i] == 1 for i in range(len(pop.ids))]
    info = None
    for outer in range(5000):
        targets, info = [], []
        for gr in groups:
            zs = [row_z(pop, s, s["v"], p, gr["scale"]) for s in gr["specs"]]
            ses = [se_on(s["qf"], lambda v, s=s: row_z(pop, s, v, p, gr["scale"])) if len(gr["specs"]) > 1 else 1
                   for s in gr["specs"]]
            if len(gr["specs"]) == 1:
                ses = [0 if gr["specs"][0]["qf"] is None else 1]
            zstar, inf = pool(zs, ses, method)
            inf.update(z=zs, se=ses, estimate=zstar)
            info.append(inf)
            s1 = gr["specs"][int(np.argmax(np.array(ses) == 0))] if inf.get("exact") else gr["specs"][0]
            targets.append(z_target(pop, s1, zstar, p, gr["scale"]))
        for sweep in range(10):
            for m, t in zip(masks + margins, targets + list(pop.P)):
                s = p[m].sum()
                p[m] *= t / s
                p[~m] *= (1 - t) / (1 - s)
        err = max(abs(p[m].sum() - t) for m, t in zip(masks + margins, targets + list(pop.P)))
        new = []
        for gr, inf in zip(groups, info):
            s1 = gr["specs"][0]
            new.append(z_target(pop, s1 if not inf.get("exact") else
                                gr["specs"][int(np.argmax(np.array(inf["se"]) == 0))], inf["estimate"], p, gr["scale"]))
        if err < tol and max(abs(a - b) for a, b in zip(new, targets)) < tol:
            break
    return p, groups, info


def five_population():
    d = read("diseases.csv")
    ids = [r["disease"] for r in d]
    P = [1 - np.exp(-float(r["value"])) if r["measure"] == "incidence_rate" else float(r["value"]) for r in d]
    pop = Pop(ids, P)
    specs = [spec(pop, r) for r in read("associations.csv")]
    return pop, specs


# ---------------------------------------------------------------- impacts
def impact_spec(pop, r):
    g = split(r["disease"])
    apos, aneg = pop.among(r.get("among", ""), g)
    adj = split(r.get("adjusted_for", ""))
    return dict(row=r["row"], g=g, term=";".join([d for d in pop.ids if d in g]), apos=apos, aneg=aneg,
                all=adj == ["all"], adj=[] if adj == ["all"] else [x for x in adj if x in pop.ids],
                m=(r.get("measure") or "").strip(), v=float(r["value"]), qf=quantile_fn(r))


def terms_of(pop, specs):
    combos = []
    for s in specs:
        if len(s["g"]) > 1 and s["term"] not in combos:
            combos.append(s["term"])
    return pop.ids + combos


def additive_coefs(pop, p, s, terms):
    F = np.column_stack([pop.allm(split(t)) for t in terms]).astype(float)
    A = pop.allm(s["apos"]) & pop.nonem(s["aneg"])
    XG = pop.allm(s["g"])
    if not s["all"] and not s["adj"]:
        w1, w0 = p * (A & XG), p * (A & pop.nonem(s["g"]))
        return (w1 @ F) / w1.sum() - (w0 @ F) / w0.sum()
    cols = [XG[A]]
    names = [s["term"]]
    if len(s["g"]) > 1:
        for d in s["g"]:
            cols.append(pop.col(d)[A]); names.append(d)
        for k in range(2, len(s["g"])):
            for sub in itertools.combinations(s["g"], k):
                cols.append(pop.allm(list(sub))[A]); names.append(";".join(sub))
    others = ([d for d in pop.ids if d not in s["g"] + s["apos"] + s["aneg"]] if s["all"] else s["adj"])
    for d in others:
        if d not in names:
            cols.append(pop.col(d)[A]); names.append(d)
    if s["all"]:
        for t in terms:
            if t not in names:
                cols.append(pop.allm(split(t))[A]); names.append(t)
    Z = np.column_stack([np.ones(A.sum())] + [c.astype(float) for c in cols])
    sw = np.sqrt(p[A])
    B = np.linalg.lstsq(Z * sw[:, None], F[A] * sw[:, None], rcond=None)[0]
    return B[1]


def wls(Am, y, se):
    """Weighted least squares with exact rows (se = 0) as constraints."""
    se = np.array(se, float)
    ex = se == 0
    if Am.shape[0] == Am.shape[1]:
        return np.linalg.solve(Am, y)
    AE, yE = Am[ex], y[ex]
    # Null-space method.
    U, d, Vt = np.linalg.svd(AE)
    r = int((d > 1e-9 * d.max()).sum())
    th0 = Vt[:r].T @ ((U[:, :r].T @ yE) / d[:r])
    N = Vt[r:].T
    AU, yU, w = Am[~ex], y[~ex], 1 / se[~ex]
    phi = np.linalg.lstsq((AU @ N) * w[:, None], (yU - AU @ th0) * w, rcond=None)[0]
    return th0 + N @ phi


def adjust_additive(pop, p, rows):
    specs = [impact_spec(pop, r) for r in rows]
    terms = terms_of(pop, specs)
    Am = np.array([additive_coefs(pop, p, s, terms) for s in specs])
    y = np.array([s["v"] for s in specs])
    se = [se_on(s["qf"], lambda v: v) for s in specs]
    th = wls(Am, y, se)
    Pm = pop.cells.T @ p
    pT = np.array([p[pop.allm(split(t))].sum() for t in terms])
    contrib = {d: Pm[i] * th[i] for i, d in enumerate(pop.ids)}
    for j, t in enumerate(terms[len(pop.ids):], start=len(pop.ids)):
        mem = split(t)
        for d in mem:
            contrib[d] += th[j] * pT[j] / len(mem)
    resid = (y - Am @ th) / np.where(np.array(se) > 0, se, np.inf)
    return dict(zip(terms, th)), contrib, sum(contrib.values()), resid, Am, y, se


# ----------------------------------------------------------------- events
def event_setup(pop, p, s, terms):
    if s["all"] and s["m"] in ("HR", "rate_ratio"):
        return dict(own=terms.index(s["term"]))
    A = pop.allm(s["apos"]) & pop.nonem(s["aneg"]) & (p > 0)
    k = len(s["g"])
    sub = sum((pop.col(d).astype(int) << j) for j, d in enumerate(s["g"]))
    full = 2 ** k - 1
    adjusted = s["all"] or bool(s["adj"])
    if adjusted:
        S = [d for d in pop.ids if d not in s["g"] + s["apos"] + s["aneg"]] if s["all"] else s["adj"]
        strat = sum((pop.col(d).astype(int) << j) for j, d in enumerate(S)) if S else np.zeros(len(p), int)
        need = list(range(full + 1))
    else:
        strat = np.zeros(len(p), int)
        need = [0, full]
    groups = {}
    for st in np.unique(strat[A]):
        cells = {u: np.where(A & (strat == st) & (sub == u))[0] for u in need}
        if all(p[c].sum() > 0 for c in cells.values()):
            groups[st] = cells
    return dict(own=None, groups=groups, k=k, full=full, adjusted=adjusted)


def event_value(su, m, theta, mult, R, p):
    if su["own"] is not None:
        return theta[su["own"]]
    val = mult if m in ("HR", "rate_ratio") else R
    num_, den = 0.0, 0.0
    for st, cells in su["groups"].items():
        d = {u: p[c].sum() for u, c in cells.items()}
        e = {}
        for u, c in cells.items():
            mm = (p[c] * val[c]).sum() / d[u]
            e[u] = {"HR": np.log, "rate_ratio": np.log, "RR": np.log,
                    "OR": lambda x: np.log(x / (1 - x)), "RD": lambda x: x}[m](mm)
        if su["adjusted"]:
            contr = sum((-1) ** (su["k"] - bin(u).count("1")) * e[u] for u in e)
        else:
            contr = e[su["full"]] - e[0]
        w = 1 / sum(1 / x for x in d.values())
        num_ += w * contr
        den += w
    return num_ / den


def baseline_h0(mult, p, r):
    return brentq(lambda h: (p * (1 - np.exp(-h * mult))).sum() - r, 0, 50 / mult.min(), xtol=1e-16)


def adjust_event(pop, p, rows, risk):
    specs = [impact_spec(pop, r) for r in rows]
    terms = terms_of(pop, specs)
    F = np.column_stack([pop.allm(split(t)) for t in terms]).astype(float)
    sus = [event_setup(pop, p, s, terms) for s in specs]
    target = np.array([s["v"] if s["m"] == "RD" else np.log(s["v"]) for s in specs])

    def resid(theta):
        mult = np.exp(F @ theta)
        h0 = baseline_h0(mult, p, risk)
        R = 1 - np.exp(-h0 * mult)
        return np.array([event_value(su, s["m"], theta, mult, R, p) for su, s in zip(sus, specs)]) - target

    sol = least_squares(resid, np.zeros(len(terms)), xtol=1e-15, ftol=1e-15, gtol=1e-15)
    theta = sol.x
    mult = np.exp(F @ theta)
    h0 = baseline_h0(mult, p, risk)
    r0 = 1 - np.exp(-h0)
    # Shapley allocation of the attributable risk over cells.
    phi = {d: 0.0 for d in pop.ids}
    for c, pc in zip(pop.cells, p):
        present = [d for d in pop.ids if c[pop.ids.index(d)] == 1]
        if not present:
            continue

        def loss(subset):
            lp = sum(theta[j] for j, t in enumerate(terms) if set(split(t)) <= set(subset))
            return 1 - np.exp(-h0 * np.exp(lp)) - r0

        k = len(present)
        for d in present:
            rest = [x for x in present if x != d]
            for size in range(k):
                for sub in itertools.combinations(rest, size):
                    w = np.math.factorial(size) * np.math.factorial(k - size - 1) / np.math.factorial(k) \
                        if hasattr(np, "math") else None
                    import math
                    w = math.factorial(size) * math.factorial(k - size - 1) / math.factorial(k)
                    phi[d] += pc * w * (loss(list(sub) + [d]) - loss(list(sub)))
    return dict(zip(terms, np.exp(theta))), risk - r0, phi, np.abs(sol.fun).max()


# ------------------------------------------------------------------- run
if __name__ == "__main__":
    pop, specs = five_population()
    p, groups, info = fit(pop, specs)
    print("Five diseases (inst/extdata/five_diseases)")
    for gr, inf in zip(groups, info):
        if len(gr["specs"]) > 1:
            print(f"  pooled LAM-MAS odds ratio: {np.exp(inf['estimate']):.10f}  Q {inf['Q']:.10f} tau2 {inf['tau2']:.10f}")
    LMS = pop.allm(["LAM", "MAS", "SCK"])
    print(f"  P(LAM, MAS, SCK) = {p[LMS].sum():.12f}")
    # OR ratio reproduced
    s12 = [s for s in specs if s["m"] == "OR_ratio"][0]

    def stratum_or(pres):
        A = pop.allm(["SCK"]) if pres else pop.nonem(["SCK"])
        pa = p[A].sum()
        return odds_ratio(p[A & pop.allm(["LAM", "MAS"])].sum() / pa, p[A & pop.col("LAM")].sum() / pa,
                          p[A & pop.col("MAS")].sum() / pa)
    print(f"  OR ratio at the fit: {stratum_or(True) / stratum_or(False):.12f}")
    Pm = pop.cells.T @ p
    for f in ("impacts_yield.csv", "impacts_calving_interval.csv", "impacts_welfare.csv"):
        b, contrib, total, resid, Am, y, se = adjust_additive(pop, p, read(f))
        print(f"  {f}: adjusted", {k: round(v, 10) for k, v in b.items()})
        print(f"     total {total:.10f}; contributions", {k: round(v, 10) for k, v in contrib.items()})
        if Am.shape[0] > Am.shape[1]:
            print(f"     std residuals {np.round(resid, 10)}; Q {np.nansum(resid[np.isfinite(resid)] ** 2):.10f}")
    for risk in (0.25,):
        hr, ar, phi, res = adjust_event(pop, p, read("culling.csv"), risk)
        print(f"  culling (overall risk {risk}): adjusted HRs", {k: round(v, 10) for k, v in hr.items()})
        print(f"     attributable {ar:.10f}; Shapley", {k: round(v, 10) for k, v in phi.items()},
              f"(sum {sum(phi.values()):.10f}; residual {res:.1e})")


def bob_three_way():
    """Section 4 of test-review-v03.R: three diseases, p = 0.5, no pair row,
    one OR_ratio row (a with b among c) of 4; crude impacts 1, 2, 3."""
    pop = Pop(["a", "b", "c"], [0.5, 0.5, 0.5])
    r = dict(row=1, disease="a", associated_with="b", among="c", adjusted_for="", measure="OR_ratio",
             value="4", dist="")
    p, groups, info = fit(pop, [spec(pop, r)])
    rows = [dict(row=i + 1, disease=d, among="", adjusted_for="", value=str(v), dist="")
            for i, (d, v) in enumerate(zip("abc", (1, 2, 3)))]
    b, contrib, total, *_ = adjust_additive(pop, p, rows)
    print("Bob's three-way (OR_ratio 4, all pairs unknown): P(a,b) =", round(p[pop.allm(["a", "b"])].sum(), 12),
          " total", round(total, 12), " adjusted", {k: round(v, 12) for k, v in b.items()})


if __name__ == "__main__":
    bob_three_way()
