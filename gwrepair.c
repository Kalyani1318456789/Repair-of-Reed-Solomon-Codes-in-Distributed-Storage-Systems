#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>


typedef struct {
    int       q;        /* Base prime                  */
    int       m;        /* Extension degree            */
    long long size;     /* Field size = q^m            */
    int      *irr;      /* Irreducible polynomial      */
} Field;


/* ============================================================================
 *  Utility: Integer Exponentiation and Primality
 * ============================================================================ */

static long long ipow(long long base, int exp)
{
    long long result = 1;
    while (exp-- > 0)
        result *= base;
    return result;
}

static int is_prime(int n)
{
    if (n < 2) return 0;
    for (int i = 2; (long long)i * i <= n; i++) {
        if (n % i == 0)
            return 0;
        return 1;
    }
    return -1;
}

/* ============================================================================
 *  Vector Encoding / Decoding  (element <-> coefficient array)
 * ============================================================================ */

static void to_vec(const Field *F, long long a, int *v)
{
    for (int i = 0; i < F->m; i++) {
        v[i] = (int)(a % F->q);
        a   /= F->q;
    }
}

static long long from_vec(const Field *F, const int *v)
{
    long long x = 0, p = 1;
    for (int i = 0; i < F->m; i++) {
        x += v[i] * p;
        p *= F->q;
    }
    return x;
}


/* ============================================================================
 *  GF(q^m) Arithmetic
 * ============================================================================ */

static long long fadd(const Field *F, long long a, long long b)
{
    int *A = malloc(F->m * sizeof(int));
    int *B = malloc(F->m * sizeof(int));
    int *R = malloc(F->m * sizeof(int));

    to_vec(F, a, A);
    to_vec(F, b, B);
    for (int i = 0; i < F->m; i++)
        R[i] = (A[i] + B[i]) % F->q;

    long long r = from_vec(F, R);
    free(A); free(B); free(R);
    return r;
}

static long long fsub(const Field *F, long long a, long long b)
{
    int *A = malloc(F->m * sizeof(int));
    int *B = malloc(F->m * sizeof(int));
    int *R = malloc(F->m * sizeof(int));

    to_vec(F, a, A);
    to_vec(F, b, B);
    for (int i = 0; i < F->m; i++)
        R[i] = ((A[i] - B[i]) % F->q + F->q) % F->q;

    long long r = from_vec(F, R);
    free(A); free(B); free(R);
    return r;
}

static long long fmul(const Field *F, long long a, long long b)
{
    if (!a || !b) return 0;

    int m = F->m, q = F->q;
    int *A = calloc(m,     sizeof(int));
    int *B = calloc(m,     sizeof(int));
    int *p = calloc(2 * m, sizeof(int));

    to_vec(F, a, A);
    to_vec(F, b, B);

    /* Polynomial multiplication */
    for (int i = 0; i < m; i++)
        for (int j = 0; j < m; j++)
            p[i + j] += A[i] * B[j];

    for (int i = 0; i < 2 * m; i++)
        p[i] = ((p[i] % q) + q) % q;

    /* Reduce modulo the irreducible polynomial */
    for (int i = 2 * m - 2; i >= m; i--) {
        int c = p[i];
        if (!c) continue;
        for (int j = 0; j <= m; j++)
            p[i - m + j] = ((p[i - m + j] - c * F->irr[j]) % q + q) % q;
    }

    long long r = from_vec(F, p);
    free(A); free(B); free(p);
    return r;
}

static long long fscalar(const Field *F, int s, long long a)
{
    s = ((s % F->q) + F->q) % F->q;
    if (!s || !a) return 0;

    int *A = malloc(F->m * sizeof(int));
    int *R = malloc(F->m * sizeof(int));

    to_vec(F, a, A);
    for (int i = 0; i < F->m; i++)
        R[i] = (s * A[i]) % F->q;

    long long r = from_vec(F, R);
    free(A); free(R);
    return r;
}

static long long fpow(const Field *F, long long a, long long e)
{
    long long r = 1;
    while (e > 0) {
        if (e & 1) r = fmul(F, r, a);
        a  = fmul(F, a, a);
        e >>= 1;
    }
    return r;
}

static long long finv(const Field *F, long long a)
{
    return fpow(F, a, F->size - 2);
}

static int ftrace(const Field *F, long long x)
{
    long long s = 0, y = x;
    for (int i = 0; i < F->m; i++) {
        s = fadd(F, s, y);
        y = fpow(F, y, F->q);
    }
    int *v = malloc(F->m * sizeof(int));
    to_vec(F, s, v);
    int t = v[0];
    free(v);
    return t;
}


/* ============================================================================
 *  Irreducible Polynomial Utilities
 * ============================================================================ */

/* Returns 1 if polynomial p[0..d] has no roots in GF(q) */
static int _no_roots(const int *p, int d, int q)
{
    for (int x = 0; x < q; x++) {
        int v = 0, pw = 1;
        for (int i = 0; i <= d; i++) {
            v  = (v + (long long)p[i] * pw % q) % q;
            pw = (long long)pw * x % q;
        }
        if (v == 0) return 0;
    }
    return 1;
}

/* Returns 1 if D[0..dd] divides P[0..dp] over GF(q) */
static int _divides(const int *P, int dp, const int *D, int dd, int q)
{
    int *r = malloc((dp + 1) * sizeof(int));
    memcpy(r, P, (dp + 1) * sizeof(int));

    for (int i = dp; i >= dd; i--) {
        int c = r[i];
        if (!c) continue;
        for (int j = 0; j <= dd; j++)
            r[i - dd + j] = ((r[i - dd + j] - (long long)c * D[j]) % q + q) % q;
    }

    int ok = 1;
    for (int i = 0; i < dd; i++)
        if (r[i]) { ok = 0; break; }

    free(r);
    return ok;
}

/* Search for a monic irreducible polynomial of degree m over GF(q) */
static int find_irr_runtime(int q, int m, int *out)
{
    long long total = ipow(q, m);
    int *poly = calloc(m + 1, sizeof(int));
    int *dp2  = calloc(m + 1, sizeof(int));

    for (long long enc = 0; enc < total; enc++) {
        poly[m] = 1;
        long long tmp = enc;
        for (int i = 0; i < m; i++) {
            poly[i] = (int)(tmp % q);
            tmp    /= q;
        }

        if (!_no_roots(poly, m, q)) continue;

        int ok = 1;
        if (m >= 4) {
            for (int d = 2; d <= m / 2 && ok; d++) {
                long long tot2 = ipow(q, d);
                for (long long e2 = 0; e2 < tot2 && ok; e2++) {
                    dp2[d] = 1;
                    long long t2 = e2;
                    for (int i = 0; i < d; i++) {
                        dp2[i] = (int)(t2 % q);
                        t2    /= q;
                    }
                    if (!_no_roots(dp2, d, q)) continue;
                    if (_divides(poly, m, dp2, d, q)) ok = 0;
                }
            }
        }

        if (ok) {
            memcpy(out, poly, (m + 1) * sizeof(int));
            free(poly); free(dp2);
            return 1;
        }
    }

    free(poly); free(dp2);
    return 0;
}

/* Load irreducible polynomial into F->irr via runtime search */
static int load_irr_poly(Field *F)
{
    return find_irr_runtime(F->q, F->m, F->irr);
}


/* ============================================================================
 *  Basis Construction
 * ============================================================================ */

/* Standard polynomial basis: b[i] = X^i */
static void make_basis(const Field *F, long long *b)
{
    long long p = 1;
    for (int i = 0; i < F->m; i++) {
        b[i] = p;
        p   *= F->q;
    }
}

/* Solve  A * x = b  over GF(p) using Gauss-Jordan elimination */
static void gauss_solve(int n, int **A, int *b, int p, int *x)
{
    int **M = malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
        M[i] = malloc((n + 1) * sizeof(int));
        for (int j = 0; j < n; j++)
            M[i][j] = ((A[i][j] % p) + p) % p;
        M[i][n] = ((b[i] % p) + p) % p;
    }

    int row = 0;
    for (int col = 0; col < n; col++) {
        int sel = -1;
        for (int r = row; r < n; r++) {
            if (M[r][col]) { sel = r; break; }
        }
        if (sel < 0) continue;
        if (sel != row) { int *t = M[row]; M[row] = M[sel]; M[sel] = t; }

        /* Compute modular inverse of pivot */
        int piv = M[row][col], inv = 1, base = piv, e = p - 2;
        while (e > 0) {
            if (e & 1) inv  = (int)((long long)inv  * base % p);
            base = (int)((long long)base * base % p);
            e  >>= 1;
        }
        for (int c = col; c <= n; c++)
            M[row][c] = (int)((long long)M[row][c] * inv % p);

        for (int r = 0; r < n; r++) {
            if (r == row || !M[r][col]) continue;
            int f = M[r][col];
            for (int c = col; c <= n; c++)
                M[r][c] = ((M[r][c] - (int)((long long)f * M[row][c] % p)) % p + p) % p;
        }
        row++;
        if (row == n) break;
    }

    for (int i = 0; i < n; i++) x[i] = 0;
    for (int r = n - 1; r >= 0; r--) {
        int lead = -1;
        for (int c = 0; c < n; c++) {
            if (M[r][c]) { lead = c; break; }
        }
        if (lead < 0) continue;
        int s = M[r][n];
        for (int c = lead + 1; c < n; c++)
            s = ((s - (int)((long long)M[r][c] * x[c] % p)) % p + p) % p;
        x[lead] = s;
    }

    for (int i = 0; i < n; i++) free(M[i]);
    free(M);
}

/* Compute the dual basis w.r.t. the trace bilinear form */
static void make_dual(const Field *F, const long long *basis, long long *dual)
{
    int m = F->m, q = F->q;
    int **G = malloc(m * sizeof(int *));

    for (int i = 0; i < m; i++) {
        G[i] = malloc(m * sizeof(int));
        for (int j = 0; j < m; j++)
            G[i][j] = ftrace(F, fmul(F, basis[i], basis[j]));
    }

    int *ej = calloc(m, sizeof(int));
    int *co = calloc(m, sizeof(int));

    for (int j = 0; j < m; j++) {
        memset(ej, 0, m * sizeof(int));
        ej[j] = 1;
        gauss_solve(m, G, ej, q, co);
        dual[j] = from_vec(F, co);
    }

    for (int i = 0; i < m; i++) free(G[i]);
    free(G); free(ej); free(co);
}


/* ============================================================================
 *  Reed-Solomon Encoding
 * ============================================================================ */

/* Evaluate polynomial msg[0..k-1] at each point in eval[0..n-1] (Horner's method) */
static void rs_encode(const Field *F, const long long *msg, int k,
                      const long long *eval, int n, long long *code)
{
    for (int i = 0; i < n; i++) {
        long long y = 0, a = eval[i];
        for (int j = k - 1; j >= 0; j--)
            y = fadd(F, fmul(F, y, a), msg[j]);
        code[i] = y;
    }
}


/* ============================================================================
 *  Printing
 * ============================================================================ */

static void print_vec(const Field *F, long long a)
{
    int *v = malloc(F->m * sizeof(int));
    to_vec(F, a, v);
    printf("[");
    for (int i = F->m - 1; i >= 0; i--) {
        printf("%d", v[i]);
        if (i > 0) printf(", ");
    }
    printf("]");
    free(v);
}


/* ============================================================================
 *  Guruswami-Wootters Repair
 * ============================================================================ */

static long long repair_gw(const Field *F,
                            const long long *basis, const long long *dual,
                            const long long *eval,  int n,
                            const long long *code,  int missing, int verbose)
{
    int       m     = F->m;
    int       q     = F->q;
    long long astar = eval[missing];

    /*
     * (Optional) Ratio subsymbols as in the GW description:
     *   ratio_j = f(alpha_j) / (alpha_j - alpha*)
     *   s_j[t]  = Tr(t * ratio_j)  for t = 1..q-1
     * Retained for exposition; not used in reconstruction below.
     */
    int **s = malloc(n * sizeof(int *));
    for (int j = 0; j < n; j++) {
        s[j] = calloc(q, sizeof(int));
        if (j == missing) continue;
        long long diff  = fsub(F, eval[j], astar);
        long long ratio = fmul(F, code[j], finv(F, diff));
        for (int t = 1; t < q; t++)
            s[j][t] = ftrace(F, fscalar(F, t, ratio));
    }

    /*
     * Reconstruction via the sum-zero identity:
     *   tr[i]  =  -sum_{j != missing}  Tr(dual[i] * f(alpha_j))
     *   f(a*)  =   sum_i  tr[i] * basis[i]
     */
    int *tr = calloc(m, sizeof(int));
    for (int j = 0; j < n; j++) {
        if (j == missing) continue;
        for (int i = 0; i < m; i++) {
            int contrib = ftrace(F, fmul(F, dual[i], code[j]));
            tr[i] = (tr[i] - contrib % q + q) % q;
        }
    }

    long long recon = 0;
    for (int i = 0; i < m; i++)
        recon = fadd(F, recon, fscalar(F, tr[i], basis[i]));

    for (int j = 0; j < n; j++) free(s[j]);
    free(s);
    free(tr);
    return recon;
}


/* ============================================================================
 *  Main
 * ============================================================================ */

int main(int argc, char *argv[])
{
    /* Seed RNG */
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    srand((unsigned)(ts.tv_nsec ^ ts.tv_sec));

    /* Parse or prompt for q and m */
    int q_val = 2, m_val = 3;
    if (argc >= 3) {
        q_val = atoi(argv[1]);
        m_val = atoi(argv[2]);
    } else {
        printf("GW RS Repair over GF(q^m)\n");
        printf("Enter q (prime, use q=2 for the single-subsymbol GW formula): ");
        fflush(stdout); scanf("%d", &q_val);
        printf("Enter m (>= 2):  ");
        fflush(stdout); scanf("%d", &m_val);
    }

    if (!is_prime(q_val)) { fprintf(stderr, "q=%d not prime.\n",    q_val); return 1; }
    if (m_val < 2)         { fprintf(stderr, "m must be >= 2.\n");           return 1; }

    /* Initialise field */
    Field F;
    F.q    = q_val;
    F.m    = m_val;
    F.size = ipow((long long)q_val, m_val);
    F.irr  = calloc(m_val + 1, sizeof(int));

    if (!load_irr_poly(&F)) {
        fprintf(stderr, "ERROR: no irreducible polynomial found.\n");
        return 1;
    }

    /* Build polynomial and dual bases */
    long long *basis = malloc(m_val * sizeof(long long));
    long long *dual  = malloc(m_val * sizeof(long long));
    make_basis(&F, basis);
    make_dual(&F, basis, dual);

    /*
     * GW construction parameters:
     *   n = q^m          (evaluate on all of GF(q^m))
     *   k in [1, q^m - q^(m-1)]   chosen uniformly at random
     *       = [1, n * (1 - 1/q)]
     */
    int n_rs  = (int)F.size;
    int k_max = (int)(F.size - ipow((long long)q_val, m_val - 1)); /* q^m - q^(m-1) */
    int k     = 1 + rand() % k_max;   /* uniform in [1, k_max] */

    long long *eval = malloc(n_rs * sizeof(long long));
    long long *msg  = malloc(k_max * sizeof(long long)); /* allocate max to allow reuse */
    long long *code = malloc(n_rs * sizeof(long long));

    for (int i = 0; i < n_rs; i++) eval[i] = i;
    for (int i = 0; i < k;    i++) msg[i]  = (long long)(rand() % F.size);

    rs_encode(&F, msg, k, eval, n_rs, code);

    int       missing  = rand() % n_rs;
    long long true_val = code[missing];

    /* Print field info */
    printf("\nGF(%d^%d)   n=%d (full field)   k=%d  (random in [1, %d])\n", q_val, m_val, n_rs, k, k_max);
    printf("Irreducible polynomial: p(x) = x^%d", m_val);
    for (int i = m_val - 1; i >= 0; i--)
        if (F.irr[i]) printf(" + %d*x^%d", F.irr[i], i);
    printf("\n");

    printf("Polynomial basis b[i] = X^i:   ");
    for (int i = 0; i < m_val; i++) {
        printf("b[%d]=", i); print_vec(&F, basis[i]); printf("  ");
    }
    printf("\nDual basis u[i]:               ");
    for (int i = 0; i < m_val; i++) {
        printf("u[%d]=", i); print_vec(&F, dual[i]); printf("  ");
    }
    printf("\n\n");

    /* Print codeword */
    printf("Full codeword  (c%d is the erased symbol):\n", missing);
    for (int i = 0; i < n_rs; i++) {
        if (i == missing) {
            printf("  c%d = [ERASED]\n", i);
        } else {
            printf("  c%d = ", i); print_vec(&F, code[i]); printf("\n");
        }
    }

    /* Repair */
    printf("\n--- Repair of c%d   (alpha* = %lld) ---\n\n", missing, eval[missing]);
    long long recon = repair_gw(&F, basis, dual, eval, n_rs, code, missing, 1);

    printf("\n  c%d (erased)  = ", missing); print_vec(&F, true_val);
    printf("\n  Recovered     = ");           print_vec(&F, recon);
    printf("\n  Result        : %s\n", (recon == true_val) ? "CORRECT" : "WRONG");

    /* Batch trials */
    printf("\nRunning 50 random repair trials...\n");
    int pass = 0, fail = 0;
    for (int trial = 0; trial < 50; trial++) {
        int k_trial = 1 + rand() % k_max;   /* fresh random k each trial */
        for (int i = 0; i < k_trial; i++)
            msg[i] = (long long)(rand() % F.size);
        rs_encode(&F, msg, k_trial, eval, n_rs, code);

        int       mi = rand() % n_rs;
        long long r  = repair_gw(&F, basis, dual, eval, n_rs, code, mi, 0);

        if (r == code[mi]) {
            pass++;
        } else {
            fail++;
            printf("  FAIL trial %d (k=%d): c%d expected=", trial, k_trial, mi);
            print_vec(&F, code[mi]);
            printf(" got=");
            print_vec(&F, r);
            printf("\n");
        }
    }
    printf("50 trials: %d pass, %d fail\n", pass, fail);

    /* Cleanup */
    free(eval); free(msg); free(code);
    free(basis); free(dual); free(F.irr);

    return (fail > 0) ? 1 : 0;
}
