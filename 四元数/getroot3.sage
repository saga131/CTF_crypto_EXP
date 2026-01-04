# sagemath10.x
"""
Quaternion cube root over GF(p)

Given Q = (w,x,y,z) in GF(p)^4, find all q such that q^3 = Q.
Run with: sage -python getroot3.py

Typical use: quaternion-based crypto / CTF algebra problems.
"""


def get_root3(p, Q):
    F = GF(p)  # 定义有限域GF(p)
    Q = (F(Q[0]), F(Q[1]), F(Q[2]), F(Q[3]))  # 将四元数转换为GF(p)中的元素

    def quaternion_mult(q1, q2):
        a1, b1, c1, d1 = q1
        a2, b2, c2, d2 = q2
        scalar = a1*a2 - b1*b2 - c1*c2 - d1*d2
        i = a1*b2 + b1*a2 + c1*d2 - d1*c2
        j = a1*c2 - b1*d2 + c1*a2 + d1*b2
        k = a1*d2 + b1*c2 - c1*b2 + d1*a2
        return (scalar, i, j, k)

    def quaternion_pow(q, n):
        result = (F(1), F(0), F(0), F(0))
        while n > 0:
            if n % 2 == 1:
                result = quaternion_mult(result, q)
            q = quaternion_mult(q, q)
            n = n // 2
        return result

    # 提取目标四元数的标量部分和向量部分
    W, X, Y, Z = Q
    S = W
    V = (X, Y, Z)
    N_V = (X**2 + Y**2 + Z**2)  # 向量部分的范数
    N_Q = (W**2 + N_V)        # 目标四元数的总范数

    solutions = []

    # Step 1: 检查 N_Q 是否为三次剩余
    try:
        cube_roots_NQ = N_Q.nth_root(3, all=True)
    except ValueError:
        cube_roots_NQ = []

    for n in cube_roots_NQ:
        # Step 2: 解三次方程 4*N_V*k**3 -3*n*k +1 = 0
        if N_V == 0:
            # 处理纯标量情况
            if X == 0 and Y == 0 and Z == 0:
                try:
                    s_roots = S.nth_root(3, all=True)
                    solutions.extend( (s, F(0), F(0), F(0)) for s in s_roots )
                except:
                    pass
            continue
        
        # R.<k> = PolynomialRing(F)
        R = PolynomialRing(F, 'k')
        k = R.gen()
        eq = 4*N_V * k**3 - 3*n * k + 1
        k_candidates = eq.roots(multiplicities=False)
        
        for k in k_candidates:
            # Step 3: 计算 s² = n - k²*N_V
            s_sq = n - k**2 * N_V
            if s_sq == 0:
                s_candidates = [F(0)]
            else:
                if not s_sq.is_square():
                    continue
                s_candidates = s_sq.sqrt(all=True)
            
            for s in s_candidates:
                # Step 4: 验证标量方程 s*(n -4k²*N_V) ≡ S mod p
                lhs = s * (n - 4*k**2*N_V)
                if lhs == S:
                    q = (s, k*X, k*Y, k*Z)
                    # 验证 q**3 是否等于 Q（避免计算误差）
                    if quaternion_pow(q, 3) == Q:
                        solutions.append(q)

    # 去重并输出
    solutions = list(set(solutions))
    # print(f"解为：{solutions}")
    return solutions

if __name__ == "__main__":
    p = 63173373914948586508761871207488662566773264479285518327131522282352053209317
    Q = (36698564177888078258192095739455152652959860052111216061091759447957860686074, 17870807869940687395361314550407377371850625515573380948432760072344080142389, 28335490245070169116781105091378482201161610915164600915589821149813685522901, 11951863920094324549214074577482301476865489472163720590246328864154628320061)
    sols = get_root3(p, Q)
    print(f"[+] solutions: {len(sols)}")
    for q in sols:
        print(q)