#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

long long evaluate(long long m[], int k, long long x, long long p) {
    long long result = 0;
    long long x_power = 1;

    for (int i = 0; i < k; i++) {
        result = (result + m[i] * x_power) % p;
        x_power = (x_power * x) % p;
    }

    return result;
}

void create_random_erasures(long long c[], int n, int k) {
    
    int max_erasures = n - k;

    if (max_erasures <= 0)
        return ;

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

long long modInverse(long long a, long long p) {
    return power(a, p - 2, p);
}

long long lagrange_evaluate(long long x[], long long y[], int k,
                             long long X, long long p) {

    long long PX = 0;

    for (int i = 0; i < k; i++) {

        long long Li = 1;

        for (int j = 0; j < k; j++) {
            if (j != i) {

                long long numerator = (X - x[j] + p) % p;
                long long denominator = (x[i] - x[j] + p) % p;
                long long inv = modInverse(denominator, p);

                Li = (Li * numerator) % p;
                Li = (Li * inv) % p;
            }
        }

        PX = (PX + y[i] * Li) % p;
    }

    return PX;
}

void decode_and_recover(long long C[], int n, int k) {

    long long p = n + 1;

    long long x[k], y[k];
    int count = 0;

    // Collect first k non-erased points
    for (int i = 0; i < n && count < k; i++) {
        if (C[i] != -1) {
            x[count] = i + 1;   // evaluation points were 1..n
            y[count] = C[i];
            count++;
        }
    }

    if (count < k) {
        printf("Not enough points to decode.\n");
        return;
    }

    // Final polynomial coefficients (degree k-1)
    long long poly[k];
    for (int i = 0; i < k; i++)
        poly[i] = 0;

    // Build polynomial using Lagrange expansion
    for (int i = 0; i < k; i++) {

        // Start numerator polynomial = 1
        long long temp[k];
        for (int t = 0; t < k; t++) temp[t] = 0;
        temp[0] = 1;

        long long denom = 1;

        for (int j = 0; j < k; j++) {

            if (j == i) continue;

            // Multiply temp by (x - x[j])
            long long newtemp[k];
            for (int t = 0; t < k; t++) newtemp[t] = 0;

            for (int d = 0; d < k-1; d++) {
                if (temp[d] != 0) {
                    newtemp[d+1] = (newtemp[d+1] + temp[d]) % p;
                    newtemp[d] = (newtemp[d] - temp[d] * x[j] % p + p) % p;
                }
            }

            for (int t = 0; t < k; t++)
                temp[t] = newtemp[t];

            denom = (denom * (x[i] - x[j] + p) % p) % p;
        }

        long long inv = power(denom, p-2, p);
        long long factor = (y[i] * inv) % p;

        // Add scaled temp into final poly
        for (int d = 0; d < k; d++) {
            poly[d] = (poly[d] + factor * temp[d]) % p;
        }
    }

    // Print recovered message coefficients
    printf("\nRecovered Message M = [");
    for (int i = 0; i < k; i++) {
        printf("%lld", poly[i]);
        if (i != k-1) printf(",");
    }
    printf("]\n");

    // Recompute full codeword
    printf("Recovered Codeword = [");
    for (int i = 1; i <= n; i++) {

        long long val = 0;
        long long xp = 1;

        for (int d = 0; d < k; d++) {
            val = (val + poly[d] * xp) % p;
            xp = (xp * i) % p;
        }

        printf("%lld", val);
        if (i != n) printf(",");
    }
    printf("]\n");
}



int main() {
    int n, k;

    printf("Enter n and k: ");
    scanf("%d %d", &n, &k);

    long long p = n + 1;   // Field GF(n+1)

    printf("Field = GF(%lld)\n", p);
    printf("(Ensure %lld is prime)\n\n", p);

    if (k > n) {
        printf("k must be <= n\n");
        return 0;
    }

    srand(time(NULL));

    long long m[k];

    // Generate random message coefficients
    for (int i = 0; i < k; i++)
        m[i] = rand() % p;

    long long c[n];

    // Compute codeword c = [f(1), ..., f(n)]
    for (int i = 1; i <= n; i++)
        c[i - 1] = evaluate(m, k, i, p);

    // Print message matrix
    printf("m = [");
    for (int i = 0; i < k; i++) {
        printf("%lld", m[i]);
        if (i != k - 1) printf(",");
    }
    printf("]\n");

    // Print codeword matrix
    printf("c = [");
    for (int i = 0; i < n; i++) {
        printf("%lld", c[i]);
        if (i != n - 1) printf(",");
    }
    printf("]\n");

    printf("\n");
    printf("codeword subjected to erasures :");
    create_random_erasures(c,n,k);

    // Print codeword matrix sunjected to erasures 
    printf("c = [");
    for (int i = 0; i < n; i++) {
        printf("%lld", c[i]);
        if (i != n - 1) printf(",");
    }
    printf("]\n");

    decode_and_recover(c,n,k);

    return 0;
}

