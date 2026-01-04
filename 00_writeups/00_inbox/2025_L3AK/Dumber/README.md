# Dumber  - **Year**: 2025 - **Event**: L3AK - **Tags**:   异常椭圆曲线

### Problem (题面/关键参数/附件说明)   

“Don’t try to outsmart me buddy.”

task.py

```
from Crypto.Util.number import  bytes_to_long, long_to_bytes
from sage.all import *

a,b,p = ?,?,?

pt1="L3AK{test_"
pt2="flag}"

E = EllipticCurve(Zmod(p), [a, b])
p,q=E.random_element(),E.random_element()
u=bytes_to_long(pt1.encode())*p
v=bytes_to_long(pt2.encode())*q

# I will help u <3
print(p,u,q,v)
```

output.txt

```text
(103905521866731574234430443362297034336 : 116589269353056499566212456950780999584 : 1) (171660318017081135625337806416866746485 : 122407097490400018041253306369079974706 : 1) (161940138185633513360673631821653803879 : 167867902631659599239485617419980253311 : 1) (95406403280474692216804281695624776780 : 109560844064302254814641159241201048462 : 1)
```

### Idea - 入口特征： - 关键等式/结构： - 为什么可行：  

### Exploit ```bash python exp.py ```  

来源：[Dumber - L3AK2025 | Kesero](https://k3sero.github.io/posts/Dumber-L3AK2025/)

```python
from Crypto.Util.number import bytes_to_long, long_to_bytes
import math
from sage.all import *

x1 = 103905521866731574234430443362297034336
y1 = 116589269353056499566212456950780999584
x2 = 171660318017081135625337806416866746485
y2 = 122407097490400018041253306369079974706
x3 = 161940138185633513360673631821653803879
y3 = 167867902631659599239485617419980253311
x4 = 95406403280474692216804281695624776780
y4 = 109560844064302254814641159241201048462

# Verify the points of the curve.
def verify_point(x, y, a, b, p):
    left = pow(y, 2, p)
    right = (pow(x, 3, p) + a * x + b) % p
    return left == right

# Convert a field element to a p-adic number.
def _gf_to_qq(n, qq, x):
    return ZZ(x) if n == 1 else qq(list(map(int, x.polynomial())))

# Lift a point to the p-adic numbers.
def _lift(E, p, Px, Py):
    for P in E.lift_x(Px, all=True):
        if (P.xy()[1] % p) == Py:
            return P

def attack(G, P):

    E = G.curve()
    assert E.trace_of_frobenius() == 1, f"Curve should have trace of Frobenius = 1."

    F = E.base_ring()
    p = F.characteristic()
    q = F.order()
    n = F.degree()
    qq = Qq(q, names="g")

    # Section 6.1: case where n == 1

    E = EllipticCurve(qq, [_gf_to_qq(n, qq, a) + q * ZZ.random_element(1, q) for a in E.a_invariants()])
    Gx, Gy = _gf_to_qq(n, qq, G.xy()[0]), _gf_to_qq(n, qq, G.xy()[1])
    Gx, Gy = (q * _lift(E, p, Gx, Gy)).xy()
    Px, Py = _gf_to_qq(n, qq, P.xy()[0]), _gf_to_qq(n, qq, P.xy()[1])
    Px, Py = (q * _lift(E, p, Px, Py)).xy()
    l = ZZ(((Px / Py) / (Gx / Gy)) % p)

    if n > 1:
        # Section 6.2: case where n > 1
        G0 = p ** (n - 1) * G
        G0x, G0y = _gf_to_qq(n, qq, G0.xy()[0]), _gf_to_qq(n, qq, G0.xy()[1])
        G0x, G0y = (q * _lift(E, p, G0x, G0y)).xy()
        for i in range(1, n):
            Pi = p ** (n - i - 1) * (P - l * G)
            if Pi.is_zero():
                continue

            Pix, Piy = _gf_to_qq(n, qq, Pi.xy()[0]), _gf_to_qq(n, qq, Pi.xy()[1])
            Pix, Piy = (q * _lift(E, p, Pix, Piy)).xy()
            l += p ** i * ZZ(((Pix / Piy) / (G0x / G0y)) % p)

    return int(l)

k1 = bytes_to_long(b"L3AK{test_")  
k2 = bytes_to_long(b"flag}")       

# Compute terms for the curve equation
term1 = y1**2 - x1**3
term2 = y2**2 - x2**3
term3 = y3**2 - x3**3
term4 = y4**2 - x4**3

# Differences of the terms
A1 = term1 - term2
A2 = term1 - term3
A3 = term1 - term4
A4 = term3 - term4

# Differences in x-coordinates
dx12 = x1 - x2
dx13 = x1 - x3
dx14 = x1 - x4
dx34 = x3 - x4

# Expressions that are multiples of p
d1 = A1 * dx13 - A2 * dx12
d2 = A1 * dx14 - A3 * dx12
d3 = A1 * dx34 - A4 * dx12
d4 = A2 * dx14 - A3 * dx13
d5 = A2 * dx34 - A4 * dx13

# Compute GCD of the absolute values
candidate = abs(d1)
candidate = math.gcd(candidate, abs(d2))
candidate = math.gcd(candidate, abs(d3))
candidate = math.gcd(candidate, abs(d4))
candidate = math.gcd(candidate, abs(d5))

temp = candidate
for f in range(2, 1000000):
    while temp % f == 0:
        temp //= f

p = temp

num_a = ( (y1**2 - y2**2) - (x1**3 - x2**3) ) % p
denom_a = (x1 - x2) % p
inv_denom = pow(denom_a, -1, p)
a = (num_a * inv_denom) % p

b = (y1**2 - x1**3 - a * x1) % p

assert verify_point(x1, y1, a, b, p)
assert verify_point(x2, y2, a, b, p)
assert verify_point(x3, y3, a, b, p)
assert verify_point(x4, y4, a, b, p)

print(f"p = {p}")
print(f"a = {a}")
print(f"b = {b}")

E = EllipticCurve(GF(p), [a, b])

G1 = E(x1, y1)
P1 = E(x2, y2)
G2 = E(x3, y3)
P2 = E(x4, y4)

f1 = attack(G1, P1)
f2 = attack(G2, P2)

print(f"[+] Flag: {long_to_bytes(f1) + long_to_bytes(f2)}")
```



### Notes - 坑点： - 可复用板子：
