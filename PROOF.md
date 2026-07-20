# Exact derivation

## 1. Printed cyclic presentation

The symbol $X$ denotes a principally polarized abelian fourfold. We set
$x=[\Theta]$ and use the normalization $\int_Xx^4=24$.

For a transverse intersection of $k$ theta translates,

$$
\mathrm{ch}(O_k)=(1-e^{-x})^k.
$$

We apply inclusion-exclusion to $n$ cyclically indexed surfaces $Z_i=D_i\cap D_{i+1}$. The result is

$$
\mathrm{ch}(O_Z)
=n\mathrm{ch}(O_2)-n\mathrm{ch}(O_3)
-\frac{n(n-5)}2\mathrm{ch}(O_4).
$$

A partial normalization of $r$ isolated ordinary nodes adds $r[pt]=rx^4/24$
to $\mathrm{ch}(\nu_*O_{\widetilde Z})$.

## 2. The final twist variable

We keep the final theta twist as the variable $t$. We define

$$
K_{n,r,t}=
[O_X\longrightarrow\nu_*O_{\widetilde Z}]
\otimes O_X(t\Theta).
$$

We expand this expression exactly through degree four. This gives

$$
\begin{aligned}
\mathrm{ch}(K_{n,r,t})={}&1+tx
+\left(\frac{t^2}{2}-n\right)x^2\\
&+\left(\frac{t^3}{6}-nt+2n\right)x^3\\
&+\left(
\frac{t^4}{24}-\frac{nt^2}{2}+2nt
+\frac{n^2}{2}-\frac{55n}{12}-\frac r{24}
\right)x^4.
\end{aligned}
$$

## 3. Character of a point in the secant plane

We define

$$
\alpha=1-\frac d2x^2+\frac{d^2}{24}x^4,
\qquad
\beta=x-\frac d6x^3.
$$

The declared secant plane contains one rank-one class with first Chern
coefficient $t$. This class is

$$
\alpha+t\beta
=1+tx-\frac d2x^2-\frac{dt}{6}x^3
+\frac{d^2}{24}x^4.
$$

## 4. Uniqueness within integral theta twists

We match the degree-two terms. This gives

$$
n=\frac{d+t^2}{2}.
$$

We substitute this value into the degree-three equation. We move all terms to
one side. The result is

$$
(t-3)(t^2+d)=0.
$$

Because $d>0$, the factor $t^2+d$ has no real root. Therefore, $t=3$ is the unique real solution.

It is also the unique integral line-bundle twist.

We substitute $t=3$ into the degree-two result. This gives

$$
n=\frac{d+9}{2}.
$$

The degree-four equation determines $r$ uniquely:

$$
\begin{aligned}
r
&=3^4-12n\cdot3^2+48n\cdot3+12n^2-110n-d^2\\
&=81-74n+12n^2-d^2\\
&=(d+9)(2d-1).
\end{aligned}
$$

This result is the partial-normalization length printed in Section 12.

Together, the component count and node count agree with the proposed $O_X(3\Theta)$ correction.

Under this presentation and degree convention, they do not agree with the displayed $O_X(\Theta)$.

## 5. The printed and corrected characters

We insert the printed values of $n$ and $r$, and we keep $t=1$. This gives

$$
P_d=1+x-\frac{d+8}{2}x^2
+\frac{3d+28}{6}x^3
+\frac{d^2-80}{24}x^4.
$$

We tensor this class by $O_X(2\Theta)$. The result is

$$
e^{2x}P_d
=1+3x-\frac d2x^2-\frac d2x^3+\frac{d^2}{24}x^4
=\alpha+3\beta.
$$

## 6. Exact-target distinction

If the target is the older $\alpha+\beta$ character, the necessary additive correction is

$$
(\alpha+\beta)-P_d
=4x^2-\frac{2d+14}{3}x^3+\frac{10}{3}x^4.
$$

An exact derived autoequivalence of $D^b(X)$ does not produce this correction.

On an abelian fourfold, $\mathrm{td}(X)=1$. Therefore,

$$
\chi(E,E)=\int_X\mathrm{ch}(E)^\vee\mathrm{ch}(E).
$$

The self-Euler pairing of $\alpha+3\beta$ is $8d(d+9)$. The self-Euler pairing
of $\alpha+\beta$ is $8d(d+1)$.

## 7. Hodge-contraction consequence

For $c=\sum_{k=0}^4c_kx^k$, the relevant $28\times28$ contraction map consists of six copies of

$$
D(c)=
\begin{pmatrix}
c_0&2c_1&-2c_2\\
c_1&4c_2&-6c_3\\
-2c_2&-12c_3&24c_4
\end{pmatrix}.
$$

Both the printed and corrected classes have total rank $12$ and kernel dimension $16$.

The six mixed kernel witnesses change from

$$
C_{ij}-2(B_{ij}-B_{ji})-(d+4)A_{ij}
$$

to

$$
C_{ij}-dA_{ij}.
$$

The second expression is the expected secant-plane witness.

These results describe only the Hodge-side contraction.

The Section 12 theta-cycle still requires its own
$HT^2\to\mathrm{Ext}^2$ Atiyah/evaluation computation.

The geometrically different $W_2$ construction does not supply this map.
