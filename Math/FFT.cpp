/************************************************
* Author        :RunningTime
* Created Time  :2017/5/9 21:34:07
* File Name     :UVA_12298.cpp
 ************************************************/

#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef pair<int,int> pii;
const int N = 5e4 + 5;
const int INF = 0x3f3f3f3f;
const int MOD = 1e9 + 7;
const double PI = acos(-1.0);
typedef complex<double> CD;

inline void FFT(vector<CD> &a, bool inverse) {
    int n = a.size();
    for (int i=0, j=0; i<n; ++i) {
        if (j > i) swap(a[i], a[j]);
        int k = n;
        while (j & (k >>= 1)) j &= ~k;
        j |= k;
    }

    double pi = inverse ? -PI : PI;
    for (int s=1; s<n; s<<=1) {
        double alpha = pi / s;
        for (int k=0; k<s; ++k) {
            CD omegak = exp(CD(0, alpha*k));
            for (int Ek=k; Ek<n; Ek+=s<<1) {
                int Ok = Ek + s;
                CD t = omegak * a[Ok];
                a[Ok] = a[Ek] - t;
                a[Ek] += t;
            }
        }
    }
    if (inverse) for (int i=0; i<n; ++i) a[i] /= n;
}

inline vector<double> operator * (const vector<double> &v1, const vector<double> &v2) {
    int s1 = v1.size(), s2 = v2.size(), S = 2;
    while (S < s1+s2) S <<= 1;
    vector<CD> a(S, 0), b(S, 0);
    for (int i=0; i<s1; ++i) a[i] = v1[i];
    FFT(a, false);
    for (int i=0; i<s2; ++i) b[i] = v2[i];
    FFT(b, false);
    for (int i=0; i<S; ++i) a[i] *= b[i];
    FFT(a, true);
    vector<double> ret(s1+s2-1);
    for (int i=0; i<s1+s2-1; ++i) ret[i] = a[i].real();
    return ret;
}

bool is_prime[N];
void prime_table(int n) {
    memset(is_prime, true, sizeof(is_prime));
    is_prime[0] = is_prime[1] = false;
    int m = (int)sqrt(n+0.5);
    for (int i=2; i<=m; ++i) if (is_prime[i]) {
        for (int j=i*i; j<=n; j+=i) is_prime[j] = false;
    }
}

const char* pokers = "SHCD";
int idx(char p) {
    return strchr(pokers, p) - pokers;
}

bool lost[4][N];

int main() {
    prime_table(50000);
    int a, b, c;
    while (scanf("%d%d%d", &a, &b, &c) == 3 && a) {
        memset(lost, false, sizeof(lost));
        for (int i=0; i<c; ++i) {
            int d; char s[5];
            scanf("%d%s", &d, &s);
            lost[idx(s[0])][d] = true;
        }
        vector<double> ans(1, 1), poly;
        for (int i=0; i<4; ++i) {
            poly.clear();
            poly.resize(b+1, 0);
            for (int j=4; j<=b; ++j) {
                if (!is_prime[j] && !lost[i][j]) poly[j] = 1.0;
            }
            ans = ans * poly;
            ans.resize(b+1);
        }
        for (int i=a; i<=b; ++i) {
            printf("%.0f\n", fabs(ans[i]));
        }
        puts("");
    }
    return 0;
}
