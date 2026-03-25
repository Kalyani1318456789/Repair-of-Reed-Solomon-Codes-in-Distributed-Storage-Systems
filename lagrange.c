#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/* ================= PRIME CHECK ================= */
int is_prime(long long n) {
    if (n < 2) return 0;
    for (long long i = 2; i * i <= n; i++) {
        if (n % i == 0) return 0;
    }
    return 1;
}

/* ================= FAST POWER ================= */
long long power(long long base, long long exp, long long p) {
    long long result = 1;
    base %= p;
    while (exp > 0) {
        if (exp % 2)
            result = (result * base) % p;
        base = (base * base) % p;
        exp /= 2;
    }
    return result;
}

/* ================= MOD INVERSE ================= */
long long modInverse(long long a, long long p) {
    if (a == 0) {
        printf("Error: Division by zero!\n");
        exit(1);
    }
    return power(a, p - 2, p);  // only valid if p is prime
}

/* ================= POLY EVALUATION ================= */
long long evaluate(long long m[], int k, long long x, long long p) {
    long long result = 0, xp = 1;

    for (int i = 0; i < k; i++) {
        result = (result + m[i] * xp) % p;
        xp = (xp * x) % p;
    }

    return result;
}

/* ================= RANDOM ERASURES ================= */
void create_random_erasures(long long c[], int n, int k) {
    int max_erasures = n - k;
    if (max_erasures <= 0) return;

    int t = rand() % max_erasures + 1;

    int count = 0;
    while (count < t) {
        int index = rand() % n;
        if (c[index] != -1) {
            c[index] = -1;
            count++;
        }
    }
}

/* ================= DECODE ================= */
void decode_and_recover(long long C[], int n, int k) {

    long long p = n + 1;

    long long x[k], y[k];
    int count = 0;

    /* Collect k valid points */
    for (int i = 0; i < n && count < k; i++) {
        if (C[i] != -1) {
            x[count] = i + 1;
            y[count] = C[i];
            count++;
        }
    }

    if (count < k) {
        printf("Not enough points to decode.\n");
        return;
    }

    long long poly[k];
    for (int i = 0; i < k; i++) poly[i] = 0;

    /* Lagrange interpolation */
    for (int i = 0; i < k; i++) {

        long long temp[k];
        for (int t = 0; t < k; t++) temp[t] = 0;
        temp[0] = 1;

        long long denom = 1;

        for (int j = 0; j < k; j++) {
            if (j == i) continue;

            long long newtemp[k];
            for (int t = 0; t < k; t++) newtemp[t] = 0;

            for (int d = 0; d < k - 1; d++) {
                if (temp[d] != 0) {
                    newtemp[d + 1] = (newtemp[d + 1] + temp[d]) % p;
                    newtemp[d] = (newtemp[d] - temp[d] * x[j] % p + p) % p;
                }
            }

            for (int t = 0; t < k; t++)
                temp[t] = newtemp[t];

            long long diff = (x[i] - x[j] + p) % p;
            denom = (denom * diff) % p;
        }

        long long inv = modInverse(denom, p);
        long long factor = (y[i] * inv) % p;

        for (int d = 0; d < k; d++) {
            poly[d] = (poly[d] + factor * temp[d]) % p;
        }
    }

    /* Print message */
    printf("\nRecovered Message M = [");
    for (int i = 0; i < k; i++) {
        printf("%lld", poly[i]);
        if (i != k - 1) printf(",");
    }
    printf("]\n");

    /* Recompute codeword */
    printf("Recovered Codeword = [");
    for (int i = 1; i <= n; i++) {
        long long val = 0, xp = 1;
        for (int d = 0; d < k; d++) {
            val = (val + poly[d] * xp) % p;
            xp = (xp * i) % p;
        }
        printf("%lld", val);
        if (i != n) printf(",");
    }
    printf("]\n");
}

/* ================= MAIN ================= */
int main() {
    int n, k;

    printf("Enter n and k: ");
    scanf("%d %d", &n, &k);

    long long p = n + 1;

    /* 🚨 STRICT CHECK */
    if (!is_prime(p)) {
        printf("ERROR: n+1 = %lld is NOT prime.\n", p);
        printf("This implementation only works over GF(p) where p is prime.\n");
        return 1;
    }

    printf("Field = GF(%lld)\n\n", p);

    if (k > n) {
        printf("k must be <= n\n");
        return 1;
    }

    srand(time(NULL));

    long long m[k];

    for (int i = 0; i < k; i++)
        m[i] = rand() % p;

    long long c[n];

    for (int i = 1; i <= n; i++)
        c[i - 1] = evaluate(m, k, i, p);

    printf("Message m = [");
    for (int i = 0; i < k; i++) {
        printf("%lld", m[i]);
        if (i != k - 1) printf(",");
    }
    printf("]\n");

    printf("Codeword c = [");
    for (int i = 0; i < n; i++) {
        printf("%lld", c[i]);
        if (i != n - 1) printf(",");
    }
    printf("]\n");

    create_random_erasures(c, n, k);

<<<<<<< HEAD
    // Print codeword matrix subjected to erasures 
    printf("c = [");
=======
    printf("\nAfter erasures c = [");
>>>>>>> edb3647 (final changes)
    for (int i = 0; i < n; i++) {
        printf("%lld", c[i]);
        if (i != n - 1) printf(",");
    }
    printf("]\n");

    decode_and_recover(c, n, k);

    return 0;
}
