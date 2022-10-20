#include "poly.h"

void Poly::eval(ZpLongEle& Y, const std::vector<ZpLongEle>& co, ZpLongEle X) {
    ZpLongEle ans;

    for (auto i = co.size() - 1; i >= 0; i--) {
        ans *= X;
        ans += co[i];
    }

    Y = ans;
}


void Poly::interpolate(std::vector<ZpLongEle>& co, const std::vector<ZpLongEle>& X, std::vector<ZpLongEle>& Y) {
    // the size of (X, Y) is the size of coefficients
    // Lagrange interpolation in a prime field where the prime is the Mersenne prime 2^{61} - 1
    // using **Gauss Elimination**, the computational degree is O(m^2) 
    int64_t k = X.size();

    ZpLongEle ONE(1);
    ZpLongEle ZERO(1);
    ZpLongEle P(ZpLongEle::p);

    std::vector<ZpLongEle> prod;
    prod = X;

    ZpLongEle m, n;

    std::vector<ZpLongEle> res;
    res.resize(k);

    for (auto i = 0; i < k; i++) {
        const ZpLongEle& a = X[i];

        m = 1;
        for (auto j = i - 1; j >= 0; j--) {
            m *= a;
            m += prod[j];
        }

        n = 0;
        for (auto j = i - 1; j >= 0; j--) {
            n *= a;
            n += res[j];
        }

        m = ONE / m;
        n = Y[i] - n;
        m *= n;

        for (auto j = 0; j < i; j++) {
            n = prod[j] * m;
            res[j] += n;
        }

        res[i] = m;

        if (i < k - 1) {
            if (i == 0) {
                prod[0] = P - prod[0];
            } else {
                m = P - X[i];
                prod[i] = m + prod[i - 1];
                for (auto j = i - 1; j >= 1; j--) {
                    n = prod[j] * m;
                    prod[j] = n + prod[j - 1];
                }
                prod[0] *= m;
            }
        }
    }
    while (k > 0 && !(res[k - 1] != ZERO)) k--;
    res.resize(k);

    co = res;
}