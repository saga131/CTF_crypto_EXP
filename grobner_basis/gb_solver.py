# sagemath 10.x

from itertools import permutations
from math import gcd
import sys

# =========================
# 0. 基本参数
# =========================
n = Integer(
    51346123992014031559655175082697970684167395475333311304964644302342195696224567294295344310402829868055631397452626012477836153077911918863960702413166012975512610833278396928183510844680256458504910587654478555702162520012904675084507891264626440491149901037372986379794780581794519575982104290783154358821)
# TODO: 模数
Zn = Zmod(n)

c1 = Integer(
    24319315677447194051443187371747479255067743760096699350944397167371168114725720790247700427952741142763441610892426608864047529341657185912118909492322774662161106498801616996096683224771572769829837624851608134358304727200743886784288098480619776052412530055546018020194040544627681873207943033870007025748)
c2 = Integer(
    39786564078572144284506485299728758988286521624370259031809631236140161282082325725135140316573904105763232744807968622765903200810033097797594471787169529885312330231639058934094205131482471767673473932389152433433898410194395448708947147212928341194143100969219109312823857009215175452177730385538708060222)
c3 = Integer(
    21891565045903053617066884716022892658769628974027471361516291164611317730287631893135184181539960500062549103756684176011544313547807404396373504189128374047522742778262347183232492997780550702287363024839047277383860212235327923210729120974240292438137529506434639962960659927050629299924917847027312495548)

vars_all = ['a', 'x']  # TODO: 所有变量名
target = 'x'  # TODO: 要解的变量


# =========================
# 1. 写方程组（poly == 0）
# =========================
def build_equations(R, gens):
    """
    返回 eqs = [poly1, poly2, ...]
    """
    a = gens.get('a', None)
    x = gens.get('x', None)
    F1 = (x - 13 * a) ** 7 - 1337
    F2 = (x + 37 * a) ** 5 + 7331
    F3 = (x - 1337) ** 137 + 7331
    eqs = [
        F1 ** 3 - Zn(c1),
        F2 ** 3 - Zn(c2),
        F3 ** 3 - Zn(c3)  # TODO: 填你的多项式
    ]
    return eqs


# =========================
# 工具函数
# =========================
def is_univariate_in(p, v):
    return p.variables() == (v,)


def to_univariate(p, target_name):
    Rx = PolynomialRing(Zn, target_name)
    return Rx(p.univariate_polynomial())


def solve_linear_univariate(px1):
    """
    解一元一次多项式：u*x + v = 0
    返回：[root]（Zmod(n) 元素）
    """
    u = px1.coefficient(1)
    v = px1.coefficient(0)

    try:
        return [(-v) * u ** (-1)]
    except Exception:
        g = gcd(Integer(lift(u)), n)
        print("[!] univariate linear: u not invertible, gcd(u,n) =", g)
        if 1 < g < n:
            print("[+] non-trivial factor found:", g, n // g)
        raise


def solve_linear_multivariate(p, target_var):
    """
    解多元多项式里“关于 target_var 的一次方程”
    形如：u*target_var + h = 0
    返回： [root_expr]，root_expr 仍可能是多元表达式
    """
    # 系数：u
    u = p.coefficient({target_var: 1})

    # 剩余项：h = p - u*target_var
    h = p - u * target_var

    try:
        root = (-h) * u ** (-1)
        return [root]
    except Exception:
        g = gcd(Integer(lift(u)), n)
        print("[!] multivariate linear: u not invertible, gcd(u,n) =", g)
        if 1 < g < n:
            print("[+] non-trivial factor found:", g, n // g)
        raise


def solve_linear_auto(p, target_var):
    """
    自动判断多元 / 一元，并调用对应解法
    """
    parent = p.parent()

    # 一元多项式环
    if parent.ngens() == 1:
        return solve_linear_univariate(p)

    # 多元多项式，但 target_var 次数为 1
    if p.degree(target_var) == 1:
        return solve_linear_multivariate(p, target_var)

    raise ValueError("Polynomial is not linear in target variable")


def extract_univariate(I, target_var, elim_vars):
    # 1️⃣ 直接从 GB 找
    GB = I.groebner_basis()
    for g in GB:
        if is_univariate_in(g, target_var):
            return g, to_univariate(g, target_var)

    # 2️⃣ elimination ideal
    J = I.elimination_ideal(*elim_vars)
    for g in J.gens():
        if is_univariate_in(g, target_var):
            return g, to_univariate(g, target_var)

    return None, None


# =========================
# 2. 主逻辑：自动尝试变量顺序
# =========================
def main():
    others = [v for v in vars_all if v != target]

    print("[*] target =", target)
    print("[*] trying variable orders...")

    for perm in permutations(others):
        vars_order = list(perm) + [target]
        print("\n[*] vars order:", vars_order)

        R = PolynomialRing(Zn, vars_order, order='lex')
        gens = R.gens_dict()

        eqs = build_equations(R, gens)
        if not eqs:
            print("[!] eqs empty, fill build_equations()")
            sys.exit(1)

        I = Ideal(eqs)

        target_var = gens[target]
        elim_vars = [gens[v] for v in vars_order if v != target]

        Px, Px1 = extract_univariate(I, target_var, elim_vars)
        if Px is None:
            print("    [-] no univariate poly found")
            continue

        print("    [+] found Px:", Px)
        print("    [+] degree:", Px1.degree())

        X = Px1.parent().gen()

        if Px1.degree() == 1:
            roots = solve_linear_auto(Px1, Px1.parent().gen())

        else:
            try:
                roots = Px1.roots(multiplicities=False)
            except Exception as e:
                print("    [!] roots() failed:", e)
                continue

        print("    [+] roots (mod n):", len(roots))
        for r in roots:
            print("        root =", r, "lift =", Integer(lift(r)))

        print("\n[✓] SUCCESS with vars order:", vars_order)
        return

    print("\n[✗] All variable orders failed")


if __name__ == "__main__":
    main()

