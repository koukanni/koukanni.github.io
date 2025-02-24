---
layout: post
title: "一類積性函數求和問題 (1)"
date: 2025-02-20
author: koukanni
categories: "tricki"
---

## $\text{Quick description}$
當你想求出某個積性函數 $f(N)$ 的前綴和，且此積性函數滿足以下條件：

- 質數的取值 $f(p)$ 可以表達成數個「完全積性函數」的和

- 質數冪的取值 $f(p^k)$ 可以被快速的計算

則可以使用 $\text{Min_25 Sieve}$ 在 $O(\frac{N^{3 / 4}}{log N})$ 的時間內求得。

<!--more-->

## $\text{Notations}$

- $p_k$：第 $k$ 個質數

- $lpf(n) := $ min{$p \gt 1 \ \land \ p \ \| \ n$}
    - 特別地，$lpf(1) = 1$

- $F(n) = \sum\limits_{i=1}^{n} f(n)$

- $F_{prime} (n) = \sum\limits_{1 \leq p \leq n \land p \in Prime} f(p)$

- $Q_n = \left\lbrace m : m = \left\lfloor \frac{n}{k} \right\rfloor, \  \forall 1 \leq k \leq n \right\rbrace$


## $\text{Algorithms}$

整個演算法分成兩個部分：

1. 算出 $F_{prime} (m), \ \forall m \in Q_n$
2. 透過 $F_{prime} (n)$ 算出每個 $F(m), \ \forall m \in Q_n$

### $\text{Lucy DP}$

定義 $S(x, n) = \sum\limits_{i=2}^{n} [i \in Prime \lor x \lt lpf(i)] \cdot f(i)$
> 形式化的理解就是進行埃氏篩時，篩完 $x$ 後還未篩掉的數字的函數值總和。

我們需要先找到一個完全積性函數 $g(n)$，滿足：
- $\forall p \in Prime, \ g(p) = f(p)$

接著考慮轉移式：

- 若 $x \notin Prime \lor x^2 \gt n$：

    - $S(x, n) = S(x - 1, n)$

    - 也就是 $x$ 在埃氏篩的過程中沒有篩掉任何一個數字

- 否則：
    - $S(x, n) = S(x - 1, n) - g(x) \cdot (S(x - 1, \lfloor \frac{n}{x} \rfloor) - S(x - 1, x - 1))$

    - $S(x - 1, \lfloor \frac{n}{x} \rfloor)$ 是要計算出這一輪中被 $x$ 篩掉的函數值和，但不包括小於 $x$ 的質數

      - $\pi(x) = S(x, x)$

    - 因為 $S(x - 1, \lfloor \frac{n}{x} \rfloor)$ 裡面的數字不一定都與 $x$ 互質，所以才需要要求 $g(n)$ 為完全積性函數

從上述討論得到轉移式：

$$S(x, n) = \begin{cases}
S(x-1, n) - g(x) \cdot (S(x-1, \left\lfloor \frac{n}{x} \right\rfloor) - S(x-1, x-1)) & \text{if } x \text{ is prime} \land x^2 \leq n \\
S(x-1, n) & \text{otherwise} \end{cases}$$


> 注意到我們只會用到 $Q_n$ 處的函數值。

時間複雜度為 $O(\frac{N^{3 / 4}}{log N})$

#### $\text{Implementation}$

```cpp
/*
 * V is an array of Q_n sorted in descending order.
 */
vector<i64> V;
// V = [floor(n / i) for i in range(1, sqrtN + 1)]
// V += list(range(V[-1] - 1, 0, -1))

int idx(const i64 x) {
  return the index of x in array V
}

for (const int &p : primes) {
  const i64 pp = i64(p) * p;
  const auto sp = fp[idx(p - 1)];
  for (int i = 0; i < nv; i++) {
    if (pp > V[i]) break;
    const int j = idx(V[i] / p);
    fp[i] -= fp[j] - sp;
  }
}
```

### $\text{Min_25 Sieve}$

某種程度上，$\text{Min_25 Sieve}$ 可以看成逆向的 $\text{Lucy DP}$

同樣定義 $S(x, n) = \sum\limits_{i=2}^{n} [i \in Prime \lor x \lt lpf(i)] \cdot f(i)$

- 若 $x \notin Prime \lor x^2 \gt n$：

    - $S(x - 1, n) = S(x, n)$

    - 也就是 $x$ 在埃氏篩的過程中沒有篩掉任何一個數字

- 否則，我們可以枚舉 $x$ 的指數 $c$

    - $S(x - 1, n) = S(x, n) + \sum\limits_{1 \leq c \ \land \ x^{c+1} \leq n} f(x^c) \cdot (S(x, \lfloor \frac{n}{x^c} \rfloor) - F_{prime}(x)) + f(x^{c + 1})$

    - 對於每一個大於等於 $2$ 的 $c$ 要再加上 $f(x^{c + 1})$ 是因為 $f(x)$ 本身已經包含在 $S(x, n)$ 裡面了

    - 要減掉 $F_{prime} (x)$ 是因為我們現在是在計算 $x$ 為最小質因數的數字的函數和，但是 $S(x, n / x^c)$ 卻包含 $F_{prime}(x)$

實作上跟 $\text{Lucy DP}$ 幾乎一樣，只差在一個是從小到大跑一遍質數，一個是從大到小。

時間複雜度一樣為 $O(\frac{N^{3 / 4}}{log N})$

#### $\text{Implementation}$

```cpp
auto f_k(i64 x, i64 c) {
  return the value f(x^c), where x is Prime
}

auto f_prime(i64 x) {
  return sum of f(p), where 1 <= p <= n and p is Prime
}

for (const int &p : primes | views::reverse) {
  const i64 pp = i64(p) * p;
  for (int i = 0; i < nv; i++) {
    if (pp > V[i]) break;
    i64 pc = p;
    for (int c = 1; pc <= V[i] / p; c++, pc *= p) {
      const int k = idx(V[i] / pc);
      F[i] += f_k(p, c) * (F[k] - f_prime(p)) + f_k(p, c + 1);
    }
  }
}
```

## $\text{Example 1}$

[P5325 【模板】Min_25 筛](https://www.luogu.com.cn/problem/P5325)
> 求 $f(p) = p^k \cdot (p^k - 1)$ 的前綴和

1. $f(p) = p \cdot (p - 1) = p^2 - p$
    - $id_k(n)$ 為完全積性函數
2. $f(p^k)$ 可以透過快速冪快速的算出

先用 $\text{Lucy DP}$ 算出 $id(n), \ id_2(n)$ 在質數處的前綴和，\\
再使用 $\text{Min_25 Sieve}$ 算出答案

## $\text{Example 2}$

[ABC370G Divisible by 3 ](https://atcoder.jp/contests/abc370/tasks/abc370_g)
> 求 $\sum\limits_{1 \leq i \leq N \ \land \ \sigma(i) \bmod 3 = 0} \prod\limits\binom{e_j + M - 1}{M - 1}, \text{ where } i = \prod p_j^{e_j}$

令 $f(n) = [\sigma(n) \bmod 3 = 0] \cdot \prod\limits\binom{e_i + M - 1}{M - 1}$

考慮將 $f$ 拆成兩個函數 $g, \ h$

$$
\begin{aligned}
g(n) &= M\\
h(n) &=
\begin{cases}
0, & \text{if } \sigma(n) \bmod 3 = 0 \\
M, & \text{otherwise}
\end{cases}
\end{aligned}
$$

得 $$f(n) = g(n) - h(n)$$，且 $g, \ h$ 都是積性的。

- $g(n)$ 在質數處的前綴和即為 $\pi(n)$ 乘上 $M$，容易使用 $\text{Lucy DP}$ 求得

- $h(n)$ 則需分別計算在模 $3$ 意義下同餘 $1 \text{ and } 2$ 的質數的數量，最後再乘上 $M$

### $\text{Solution}$
- [Min_25 Sieve](https://atcoder.jp/contests/abc370/submissions/62826443)
- [Black Algorithm](https://atcoder.jp/contests/abc370/submissions/62826479)
  - 不同於 $\text{Min_25 Sieve}$ 的實作方法，複雜度為 $O(N^{1-\epsilon})$，但常數極小

## $\text{Refs}$

- [oi-wiki Min_25 篩](https://oi-wiki.org/math/number-theory/min-25/#%E6%B1%82%E8%8E%AB%E6%AF%94%E4%B9%8C%E6%96%AF%E5%87%BD%E6%95%B0%E7%9A%84%E5%89%8D%E7%BC%80%E5%92%8C)
- [Lucy_Hedgehog's Post](https://projecteuler.net/best_posts=10)
- [AtCoder ABC370G Editorial](https://atcoder.jp/contests/abc370/editorial/10906)
- [The Black Algorithm](https://baihacker.github.io/main/2020/The_prefix-sum_of_multiplicative_function_the_black_algorithm.html)