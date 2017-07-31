#include <bits/stdc++.h>
using namespace std;

// 表达式求值，逆波兰式（后缀表达式）算法
// 输入(可以有空格，支持小数，实现'+-/*%')：((1+2)*5+1)/4=
// 注意：取模一定是要整型，实现版本数字全是double，强制类型转换可能导致错误
// 转换为后缀表达式： 得到：1 2 + 5 * 1 + 4 / =
// 计算后缀表达式：得到：4.000000
struct Reverse_Polish_Expression {
    stack<char> op;
    stack<double> num;
    bool error;

    int prior(char ch) {  //运算符的优先级
        switch(ch) {
            case '+':
            case '-': return 1;
            case '*':
            case '%':
            case '/': return 2;
            default:  return 0;
        }
    }
    string get_postfix(string s) {  //中缀表达式转变后缀表达式
        while (!op.empty()) op.pop();
        op.push('#');
        string ret = "";
        int len = s.length(), i = 0;
        while (i < len) {
            if(s[i] == ' ' || s[i] == '=') { i++; continue; }
            else if (s[i] == '(') op.push(s[i++]);
            else if (s[i] == ')') {
                while (op.top() != '(') {
                    ret += op.top(); ret += ' ';
                    op.pop();
                }
                op.pop(); i++;
            } else if (s[i] == '+' || s[i] == '-' || s[i] == '*' || s[i] == '/' || s[i] == '%')  {
                while (prior(op.top()) >= prior(s[i])) {
                    ret += op.top();   ret += ' ';
                    op.pop();
                }
                op.push(s[i++]);
            } else {
                while (isdigit(s[i]) || s[i] == '.') {
                    ret += s[i++];
                }
                ret += ' ';
            }
        }
        while (op.top() != '#') {
            ret += op.top(); ret += ' ';
            op.pop();
        }
        ret += '=';
        return ret;
    }
    double calc(double a, double b, char ch) {
        if (ch == '+') return a + b;
        if (ch == '-') return a - b;
        if (ch == '*') return a * b;
        if (ch == '%') return(int)((int)a %(int)b);
        if (ch == '/') {
            if (b != 0) return a / b;
            error = true; return 0;
        }
    }
    double solve(string str) {  //计算后缀表达式
        string s = get_postfix(str);
        while (!num.empty()) num.pop();
        error = false;
        int len = s.length(), i = 0;
        while (i < len)  {
            if (s[i] == ' ' || s[i] == '=') { i++; continue; }
            else if (s[i] == '+' || s[i] == '-' || s[i] == '*' || s[i] == '/' || s[i] == '%')  {
                double a = num.top(); num.pop();
                double b = num.top(); num.pop();
                num.push(calc(b, a, s[i])); i++;
            } else {
                double x = 0;
                while (isdigit(s[i])) { x = x * 10 + s[i] - '0'; i++; }
                if (s[i] == '.') {
                    double k = 10.0, y = 0;
                    i++;
                    while (isdigit(s[i])) {
                        y += ((s[i] - '0') / k);
                        i++; k *= 10;
                    }
                    x += y;
                }
                num.push(x);
            }
        }
        return num.top();
    }
}RPN;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(0);
    string str;
    getline(cin, str);
    cout << RPN.get_postfix(str) << endl;
    cout << fixed << setprecision(6) << RPN.solve(str) << endl;
    return 0;
}
