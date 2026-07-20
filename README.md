# Theta Twist

An exact, dependency-free audit of the cyclic theta-divisor construction in
Section 12 of Eyal Markman's *Secant sheaves and Weil classes on abelian
varieties* ([arXiv:2509.23403v2](https://arxiv.org/abs/2509.23403)).

## Proposed correction

Under the degree convention used elsewhere in the paper, the Section 12
object as printed does not lie in the stated secant plane. The plane-membership
claim is restored by one source-local change:

$$
\boxed{\mathcal O_X(\Theta)\longrightarrow\mathcal O_X(3\Theta)}.
$$

Equivalently, tensor the displayed object by $\mathcal O_X(2\Theta)$.
Its corrected Chern character is

$$
1+3\Theta-\frac d2\Theta^2-\frac d2\Theta^3
 +\frac{d^2}{24}\Theta^4
=\alpha+3\beta,
$$

which lies exactly in the claimed plane

$$
\operatorname{span}\left\{
e^{\sqrt{-d}\Theta},e^{-\sqrt{-d}\Theta}
\right\}.
$$

This is the unique integral theta-twist correction within the stated cyclic
presentation. If $n$ is the number of components, $r$ the
partial-normalization length, and $t$ the final integral theta twist,
coefficient matching forces

$$
n=\frac{d+t^2}{2},\qquad (t-3)(t^2+d)=0.
$$

For $d>0$, the unique real—and therefore unique integral—solution is
$t=3$. The top-degree equation then gives

$$
n=\frac{d+9}{2},\qquad r=(d+9)(2d-1),
$$

exactly the two numerical choices already printed in Section 12. Those
otherwise unusual choices provide two strong internal consistency checks for
the proposed $3\Theta$ correction.

For $d=3$, the exact computation is

$$
(1,1,-11/2,37/6,-71/24)e^{2\Theta}
=(1,3,-3/2,-3/2,3/8)=\alpha+3\beta.
$$

## Why this is not merely a fitted correction

The repository checks the following explicitly listed transforms with rational
arithmetic: theta twists, dualization, shifts/cone reversal, translation,
degree-zero line bundles, Todd/Mukai convention, multiplication pullback,
parameter relabeling, component counts, and point-normalization counts.

The degree-two equation leaves two nondual line-twist candidates. Only the
additional $+2\Theta$ twist passes degrees three and four. Dualization gives
the corresponding $\alpha-3\beta$ solution.

The normalized characteristic class

$$
\kappa(E)=\operatorname{ch}(E)
\exp\!\left(-c_1(E)/\operatorname{rank}(E)\right)
$$

is unchanged by the repair.

The Hodge-side contraction is likewise stable in rank: both the printed and
corrected classes give rank $12$ with a $16$-dimensional kernel. What
changes is the mixed witness. For $d=3$, it becomes exactly

$$
C_{ij}-3A_{ij},
$$

the secant-plane form obtained from the corrected class. This is a
consequence check, not a substitute for the still-missing Atiyah/evaluation
map.

## Separate $\alpha+\beta$ lane

The repaired Section 12 theta-cycle has character $\alpha+3\beta$, not the
specific $\alpha+\beta$ character of the older genus-four Jacobian
construction in Example 8.2.3 of
[arXiv:2502.03415v2](https://arxiv.org/abs/2502.03415).

No exact derived autoequivalence of $D^b(X)$ can identify those two classes,
because derived autoequivalences preserve the Euler pairing and their
self-Euler pairings are

$$
\chi(\alpha+3\beta,\alpha+3\beta)=8d(d+9),\qquad
\chi(\alpha+\beta,\alpha+\beta)=8d(d+1).
$$

At $d=3$, these are $288$ and $96$. The older $W_2$ construction
already realizes $\alpha+\beta$ exactly and should remain a distinct lane.

## Reproduce

No third-party Python package is required.

```bash
PYTHONPATH=src python3 -m thetatwist --d 3
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

To rebuild the committed JSON certificate and SHA-256 manifest:

```bash
PYTHONPATH=src python3 -m thetatwist --d 3 --write-results results
```

The full coefficient derivation is in [PROOF.md](PROOF.md). The generated
machine certificate is in [results/d3_certificate.json](results/d3_certificate.json).

## Scope

This repository certifies a Chern-character identity and its immediate
Hodge-contraction consequences. It does **not** prove semiregularity, provide
the missing Atiyah/evaluation map, or prove a new case of the Hodge
conjecture. Cyclic indexing of the divisors is not treated as an effective
cyclic group action without separate geometric data.

Prepared independently by Brian Ross, a graduate student at the University
of Oklahoma. This repository is not an institutional publication and does
not attribute intent to the paper's author; it presents a proposed correction
and proves its uniqueness among integral theta twists within the stated cyclic
presentation.
