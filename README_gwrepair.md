Here’s a **clean, professional README.md** you can directly paste into GitHub. It includes the fixes we discussed and is **prof-level polished** 👇

---

# Guruswami–Wootters Reed-Solomon Repair (C Implementation)

## Overview

This project implements the **Guruswami–Wootters (GW) repair scheme** for Reed-Solomon (RS) codes over finite fields ( GF(q^m) ).

The goal is to **recover a lost codeword symbol** in a distributed storage system using **trace-based projections**, without reconstructing the entire codeword.

---

## 🚀 Key Features

* Runtime construction of finite field ( GF(q^m) )
* Automatic search for irreducible polynomial (no lookup tables)
* Full Reed-Solomon encoding (evaluation over entire field)
* Dual basis computation using Gaussian elimination
* Guruswami–Wootters repair implementation
* 50 randomized trials for verification

---

## 📌 Key Parameters

| Symbol      | Meaning                    |
| ----------- | -------------------------- |
| ( q )       | Base field (must be prime) |
| ( m )       | Extension degree           |
| ( n = q^m ) | Codeword length            |
| ( k )       | Message dimension          |
| ( p(x) )    | Irreducible polynomial     |

---

## 🧠 Mathematical Background

### Finite Field ( GF(q^m) )

Elements are polynomials of degree < m over ( GF(q) ), reduced modulo an **irreducible polynomial** ( p(x) ).

---

### Reed-Solomon Encoding

A message polynomial:

[
f(x) = a_0 + a_1 x + \dots + a_{k-1} x^{k-1}
]

is evaluated over all elements of ( GF(q^m) ):

[
codeword = (f(\alpha_0), f(\alpha_1), ..., f(\alpha_{n-1}))
]

---

### 🔥 Sum-Zero Identity

When:
[
\deg(f) < q^{m-1}
]

then:

[
\sum_{\alpha \in GF(q^m)} f(\alpha) = 0
]

This works because:

[
\sum_{x \in GF(q^m)} x^e = 0 \quad \text{for } 1 \le e \le q^m - 2
]

---

### Trace Function

[
Tr(x) = x + x^q + x^{q^2} + \dots + x^{q^{m-1}}
]

Maps:

```text
GF(q^m) → GF(q)
```

---

### Dual Basis

Dual basis ( {u_i} ) satisfies:

[
Tr(b_i \cdot u_j) = \delta_{ij}
]

This allows coordinate extraction:

[
a_i = Tr(u_i \cdot a)
]

---

## ⚙️ Repair Algorithm

To recover a missing symbol ( f(\alpha^*) ):

### Step 1: Each helper node sends

[
h_{j,i} = Tr(u_i \cdot f(\alpha_j))
]

### Step 2: Accumulate

[
tr[i] = - \sum_{j \neq *} h_{j,i}
]

### Step 3: Reconstruct

[
f(\alpha^*) = \sum_i tr[i] \cdot b_i
]

---

## 🧩 Code Structure

```
Field + Utilities
    ├── ipow, is_prime
Representation
    ├── to_vec, from_vec
Field Arithmetic
    ├── fadd, fmul, fpow, finv, ftrace
Irreducible Polynomial
    ├── find_irr_runtime
Linear Algebra
    ├── gauss_solve, make_dual
RS Encoding
    ├── rs_encode
Repair
    ├── repair_gw
Main
```

---

## 🛠️ Compilation

```bash
gcc -o gw gwrepair.c -lm
```

---

## ▶️ Running

```bash
./gw 2 3
./gw 3 2
./gw 5 2
```

---

## 📊 Sample Output

```
GF(2^3) n=8 k=3
Irreducible polynomial: x^3 + x + 1

Dual basis:
u0 = 1
u1 = x²
u2 = x

Repair result:
Original = [1,0,1]
Recovered = [1,0,1]
Result: CORRECT

50 trials: 50 pass, 0 fail
```

---

## ⚠️ Limitations

* ( q ) must be **prime**
* ( m ≥ 2 )
* ( q^m ) must fit in `long long`
* Irreducible polynomial search is **brute-force** (slow for large fields)
* Supports **single erasure only**
* No bandwidth reduction (basic GW, not optimized version)

---

## 💡 Important Insight

```text
Trace breaks field elements into small pieces.
Dual basis reconstructs them exactly.
```

---

## 📚 References

* Guruswami & Wootters (2017), *Repairing Reed-Solomon Codes*
* Lidl & Niederreiter, *Finite Fields*
* Reed & Solomon (1960)

---

