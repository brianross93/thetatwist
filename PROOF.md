# Exact derivation

## 1. Printed cyclic presentation

Let $X$ be a principally polarized abelian fourfold, write
$x=[\Theta]$, and normalize $\int_Xx^4=24$. For a transverse
intersection of $k$ theta translates,

$$
\mathrm{ch}(O_k)=(1-e^{-x})^k.
$$

For $n$ cyclically indexed surfaces $Z_i=D_i\cap D_{i+1}$,
inclusion-exclusion gives

$$
\mathrm{ch}(O_Z)
=n\mathrm{ch}(O_2)-n\mathrm{ch}(O_3)
-\frac{n(n-5)}2\mathrm{ch}(O_4).
$$

Partially normalizing $r$ isolated ordinary nodes adds $r[pt]=rx^4/24$
to $\mathrm{ch}(\nu_*O_{\widetilde Z})$.

## 2. Keep the final twist variable

Set

$$
K_{n,r,t}=
[O_X\longrightarrow\nu_*O_{\widetilde Z}]
\otimes O_X(t\Theta).
$$

Exact expansion through degree four gives

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

Write

$$
\alpha=1-\frac d2x^2+\frac{d^2}{24}x^4,
\qquad
\beta=x-\frac d6x^3.
$$

A rank-one class in the declared secant plane with first Chern coefficient
$t$ is uniquely

$$
\alpha+t\beta
=1+tx-\frac d2x^2-\frac{dt}{6}x^3
+\frac{d^2}{24}x^4.
$$

## 4. Uniqueness within integral theta twists

Matching degree two yields

$$
n=\frac{d+t^2}{2}.
$$

Substituting this into the degree-three equation and moving all terms to one
side gives

$$
(t-3)(t^2+d)=0.
$$

Because $d>0$, $t^2+d$ has no real root. Thus $t=3$ is the unique real
solution and in particular the unique integral line-bundle twist.

Now

$$
n=\frac{d+9}{2}.
$$

The degree-four equation uniquely determines

$$
\begin{aligned}
r
&=3^4-12n\cdot3^2+48n\cdot3+12n^2-110n-d^2\\
&=81-74n+12n^2-d^2\\
&=(d+9)(2d-1),
\end{aligned}
$$

which is exactly the partial-normalization length printed in Section 12.

Therefore the component count and node count are simultaneously consistent
with the proposed $O_X(3\Theta)$ correction, and inconsistent with the
displayed $O_X(\Theta)$ under this presentation and degree convention.

## 5. The printed and corrected characters

Substituting the printed $n$ and $r$ but retaining $t=1$ gives

$$
P_d=1+x-\frac{d+8}{2}x^2
+\frac{3d+28}{6}x^3
+\frac{d^2-80}{24}x^4.
$$

Tensoring by $O_X(2\Theta)$ yields

$$
e^{2x}P_d
=1+3x-\frac d2x^2-\frac d2x^3+\frac{d^2}{24}x^4
=\alpha+3\beta.
$$

## 6. Exact-target distinction

If one instead requires the older $\alpha+\beta$ character, the necessary
additive correction is

$$
(\alpha+\beta)-P_d
=4x^2-\frac{2d+14}{3}x^3+\frac{10}{3}x^4.
$$

That is not achieved by an exact derived autoequivalence of $D^b(X)$. On an
abelian fourfold,
$\mathrm{td}(X)=1$, and

$$
\chi(E,E)=\int_X\mathrm{ch}(E)^\vee\mathrm{ch}(E).
$$

The two Euler pairings are $8d(d+9)$ and $8d(d+1)$, respectively.

## 7. Hodge-contraction consequence

For $c=\sum_{k=0}^4c_kx^k$, the relevant $28\times28$ contraction map
reduces to six copies of

$$
D(c)=
\begin{pmatrix}
c_0&2c_1&-2c_2\\
c_1&4c_2&-6c_3\\
-2c_2&-12c_3&24c_4
\end{pmatrix}.
$$

Both the printed and corrected classes give total rank $12$ and kernel
dimension $16$. The six mixed kernel witnesses change from

$$
C_{ij}-2(B_{ij}-B_{ji})-(d+4)A_{ij}
$$

to

$$
C_{ij}-dA_{ij}.
$$

The latter is the expected secant-plane witness. This is only the Hodge-side
contraction. The Section 12 theta-cycle still requires its own
$HT^2\to\mathrm{Ext}^2$ Atiyah/evaluation computation; it cannot
borrow that map from the geometrically different $W_2$ construction.
