/*
    *数位DP（问在区间[l, r]范围内满足x % f(x) == 0 (f(x)为所有位数相加的和)的个数．）
    *dp[len][sum][%mod][mod]一维表示长度，一维表示数位和，一维表示对mod取模，
    *该mod是全局，即对所有mod试一次．判断符合的条件就是％mod的值为0，且数位和为mod．
*/
int DFS(int len, int sum, int val, int mod, bool limit) {
    if (len == -1) {
        return val == 0 && sum == mod;
    }
    int &now = dp[len][sum][val][mod];
    if (now != -1 && !limit) {
        return now;
    }
    int ret = 0;
    int d = limit ? digit[len] : 9;
    for (int i=0; i<=d; ++i)  {
        ret += DFS (len - 1, sum + i, (val * 10 + i) % mod, mod, limit && i == d);
    }
    if (!limit) {
        now = ret;
    }
    return ret;
}
int calc(int x) {
    int p = 0;
    while (x) {
        digit[p++] = x % 10;
        x /= 10;
    }
    int ret = 0;
    for (int i=1; i<=81; ++i) {
        ret += DFS (p - 1, 0, 0, i, true);
    }
    return ret;
}
