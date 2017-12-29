from protokol import *

import warnings
warnings.filterwarnings("ignore")

f = 50
OMEGA = 2 * sp.pi * f

def init_u1():
    from chyby import metex_m_3270d, metex_mxd_4660a

    u1a = pd.read_csv("../raw/u1a.txt", index_col=0, sep=r"\s*,\s*")
    u1a.U = u1a.U.apply(metex_mxd_4660a.acv)
    u1a.I = u1a.I.apply(metex_m_3270d.aca)
    u1a.P = u1a.P.apply(lambda v: uf(v, 0.075))

    u1d = pd.read_csv("../raw/u1d.txt", index_col=0, sep=r"\s*,\s*")
    u1d.U = u1d.U.apply(lambda v: uf(v, 0.004 * v + 0.5))
    u1d.I = u1d.I.apply(lambda v: uf(v, 0.004 * v + 0.005))
    u1d.P = u1d.P.apply(lambda v: uf(v, 0.008 * v + 0.01))

    return u1a, u1d


def init_u4():
    p = pd.read_csv("../raw/u4p.txt", sep=r"\s*,\s*")
    p.C = p.C.apply(lambda v: uf(v, 0.01 * v))
    p.U = p.U.apply(lambda v: uf(v, 0.004 * v + 0.5))
    p.I = p.I.apply(lambda v: uf(v, 0.004 * v + 0.005))
    p.P = p.P.apply(lambda v: uf(v, 0.008 * v + 0.01))

    s = pd.read_csv("../raw/u4s.txt", sep=r"\s*,\s*")
    s.C = s.C.apply(lambda v: uf(v, 0.01 * v))
    s.U = s.U.apply(lambda v: uf(v, 0.004 * v + 0.5))
    s.I = s.I.apply(lambda v: uf(v, 0.004 * v + 0.005))
    s.P = s.P.apply(lambda v: uf(v, 0.008 * v + 0.01))

    return p, s


def init_u5():
    df = pd.read_csv("../raw/u5.txt", sep=r"\s*,\s*")
    df.C = df.C.apply(lambda v: uf(v, 0.01 * v))
    df.U = df.U.apply(lambda v: uf(v, 0.004 * v + 0.5))
    df.I = df.I.apply(lambda v: uf(v, 0.004 * v + 0.005))
    df.P = df.P.apply(lambda v: uf(v, 0.008 * v + 0.01))
    return df


def add_ucinik_fi(df):
    df["ucinik"] = df.P / (df.U * df.I)
    df.ucinik = df.ucinik.apply(lambda v: v if v <=0.99999 else uf(0.99999, v.s))
    df["abs_fi"] = unumpy.degrees(unumpy.arccos(df.ucinik))

    return df


def u3_seriove(df):
    U, I, _, _, fi = df.loc["L", :]

    R = (U / I) * (1 / um.sqrt(1 + sp.tan(fi.n)**2))
    L = (U / (OMEGA * I)) * (sp.sqrt(sp.tan(fi.n)**2 / (1 + sp.tan(fi.n)**2)))
    return R, L


def u3_paralelni(df):
    U, I, _, _, fi = df.loc["L", :]

    R = (U / I) * (um.sqrt(1 + sp.tan(fi.n)**2))
    L = (U / (OMEGA * I)) * (sp.sqrt((1 + sp.tan(fi.n)**2) / sp.tan(fi.n)**2))

    return R, L


def plot_u4s(df):
    fig = plt.figure()
    ax1 = plt.axes()

    spl1 = Spline(sp.append(noms(df.C)[:], 11), sp.append(noms(df.ucinik)[:], 0.93))
    # spl1 = Spline(noms(df.C), noms(df.ucinik))
    # spl1.set_smoothing_factor(1)
    ax1.plot(*spl1.curve(1, 10), "--", c="grey")
    ax1.errorbar(noms(df.C), noms(df.ucinik), yerr=stds(df.ucinik), fmt="kx", elinewidth=1, capsize=2, label="$\\cos \\varphi$")

    ax2 = ax1.twinx()
    spl2 = Spline(sp.append(noms(df.C)[:], 11), sp.append(noms(df.abs_fi)[:], 18))
    spl2.set_smoothing_factor(100)
    ax2.plot(*spl2.curve(1, 10), "-.", c="grey")
    ax2.errorbar(noms(df.C), noms(df.abs_fi), yerr=stds(df.abs_fi), fmt="ko", elinewidth=1, ms=3, capsize=2, label="$\\varphi$")

    ax1.set_xlabel("$C [\\si{\\micro\\farad}]$")
    ax1.set_ylabel("$\\cos \\varphi$")
    ax2.set_ylabel("$\\varphi [\si{\\degree}]$")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    ax1.legend(h1+h2, l1+l2, loc="lower center")

    fig.tight_layout()
    fig.savefig("../plot/u4s.pdf")

    plt.show()

def plot_u4p(df):
    fig = plt.figure()
    ax1 = plt.axes()

    spl1 = Spline(sp.delete(noms(df.C), 6), sp.delete(noms(df.ucinik), 6))
    ax1.plot(*spl1.curve(), "--", c="grey")
    ax1.errorbar(noms(df.C), noms(df.ucinik), yerr=stds(df.ucinik), fmt="kx", elinewidth=1, capsize=2, label="$\\cos \\varphi$")

    ax2 = ax1.twinx()
    spl2 = Spline(sp.delete(noms(df.C), 6), sp.delete(noms(df.abs_fi), 6))
    spl2.set_smoothing_factor(1)
    ax2.plot(*spl2.curve(), "-.", c="grey")
    ax2.errorbar(noms(df.C), noms(df.abs_fi), yerr=stds(df.abs_fi), fmt="ko", elinewidth=1, ms=3, capsize=2, label="$\\varphi$")

    ax1.set_xlabel("$C [\\si{\\micro\\farad}]$")
    ax1.set_ylabel("$\\cos \\varphi$")
    ax2.set_ylabel("$\\varphi [\si{\\degree}]$")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    ax1.legend(h1+h2, l1+l2, loc="center right")

    fig.tight_layout()
    fig.savefig("../plot/u4p.pdf")

    plt.show()


def plot_u5_ucinik_fi(df):
    ux = sp.delete(noms(df.C), [19,20])
    uy = sp.delete(noms(df.ucinik), [19,20])

    fig = plt.figure()
    ax1 = plt.axes()

    spl1 = Spline(ux, uy)
    spl1.set_smoothing_factor(0.002)
    ax1.plot(*spl1.curve(), "--", c="grey")
    ax1.plot(ux, uy, "kx", label="$\\cos \\varphi$")

    ax2 = ax1.twinx()

    fx = sp.delete(noms(df.C), [19,20])
    fy = sp.delete(noms(df.abs_fi), [19,20])

    spl2 = Spline(fx, fy)
    spl2.set_smoothing_factor(90)
    ax2.plot(*spl2.curve(), "-.", c="grey")
    ax2.plot(fx, fy, "ko", ms=3, label="$\\varphi$")

    ax1.set_xlabel("$C [\\si{\\micro\\farad}]$")
    ax1.set_ylabel("$\\cos \\varphi$")
    ax2.set_ylabel("$\\varphi [\si{\\degree}]$")

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()

    ax1.legend(h1+h2, l1+l2, loc="lower right")

    fig.tight_layout()
    fig.savefig("../plot/u5uf.pdf")

    plt.show()


def plot_u5_p(df):
    fig = plt.figure()
    ax = plt.axes()

    spl = Spline(noms(df.C), noms(df.P))
    spl.set_smoothing_factor(0.002)
    ax.plot(*spl.curve(), "--", c="grey")
    ax.plot(noms(df.C), noms(df.P), "kx")

    ax.set_xlabel("$C [\\si{\\micro\\farad}]$")
    ax.set_ylabel("$P [\\si{\\watt}]$")

    fig.tight_layout()
    fig.savefig("../plot/u5p.pdf")

    plt.show()


def main():
    u1a, u1d = init_u1()
    u4p, u4s = init_u4()
    u5 = init_u5()

    u1a = add_ucinik_fi(u1a)
    u1d = add_ucinik_fi(u1d)

    u4p = add_ucinik_fi(u4p)
    u4s = add_ucinik_fi(u4s)

    u5 = add_ucinik_fi(u5)

    # plot_u4p(u4p)
    # plot_u4s(u4s)
    # plot_u5_ucinik_fi(u5)
    # plot_u5_p(u5)

    print("{}, {}".format(*u3_paralelni(u1d)))


if __name__ == '__main__':
    main()
