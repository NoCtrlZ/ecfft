import sys

p = 0x30644E72E131A029B85045B68181585D97816A916871CA8D3C208C16D87CFD47
a = 1
b = 5612291247948481584627780310922020304781354847659642188369727566000581075360

F = GF(p)
E = EllipticCurve(F, [a, b])
O = E.order()
n = O.p_primary_part(2)
log_n = n.log(2)
k = 1

G = E([0x1c30bfe9bd1623f6260d26c12f05dd8175ce795513e0bf38f3c0656f8ff34e2c, 0x2851ed59c90496afe1c47936c576ec08d158cae28c7b02381ee3da4379313853, 1])
R = E([0x13305a0c7fac8d4e49e9941334ba73f52b4975b0074a4c1943cc2220bc2c43d9, 0x1d252a06a64cfbef58d344b7260f9ab470a01a702609f044e85005efc204f710, 1])

assert (n * G).is_zero()
H = [R + i * G for i in range(2 ^ log_n)]
L = [h.xy()[0] for h in H]
S = [L[i] for i in range(0, n, 2)]
S_prime = [L[i] for i in range(1, n, 2)]

for i in range(log_n - 1, 0, -1):
    n = 1 << i
    nn = n // 2

    for iso in E.isogenies_prime_degree(2):
        psi = iso.x_rational_map()
        if len(set([psi(x) for x in S])) == nn:
            print("n: ", n, "isogeny: ", psi)
            break
    S = [psi(x) for x in S[:nn]]
    S_prime = [psi(x) for x in S_prime[:nn]]
    E = iso.codomain()
