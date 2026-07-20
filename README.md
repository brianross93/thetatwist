# Theta Twist

This repository gives an exact check of a cyclic theta-divisor construction.
It examines Section 12 of Eyal Markman's *Secant sheaves and Weil classes on
abelian varieties* ([arXiv:2509.23403v2](https://arxiv.org/abs/2509.23403)).
The calculation uses exact rational arithmetic.

## Proposed correction

Other sections of the paper use this degree convention. Under this convention,
the printed Section 12 object does not lie in the stated secant plane.

The proposed change is:

$$
\boxed{\mathcal O_X(\Theta)\longrightarrow\mathcal O_X(3\Theta)}.
$$

After this change, the object lies in the stated secant plane.
The tensor product with $\mathcal O_X(2\Theta)$ gives the same correction.
After this correction, the Chern character is

$$
1+3\Theta-\frac d2\Theta^2-\frac d2\Theta^3
 +\frac{d^2}{24}\Theta^4
=\alpha+3\beta.
$$

This character lies exactly in the claimed plane:

$$
\mathrm{span}\left(
e^{\sqrt{-d}\Theta},e^{-\sqrt{-d}\Theta}
\right).
$$

Within the stated cyclic presentation, this is the unique integral
theta-twist correction. The variable $n$ is the number of components. The
variable $r$ is the partial-normalization length. The variable $t$ is the final
integral theta twist.

A comparison of coefficients gives

$$
n=\frac{d+t^2}{2},\qquad (t-3)(t^2+d)=0.
$$

For $d>0$, the equation has only one real solution: $t=3$. Therefore, $t=3$
is also the only integral solution. The top-degree equation then gives

$$
n=\frac{d+9}{2},\qquad r=(d+9)(2d-1).
$$

Section 12 already uses these two numerical values. Thus, both printed values
agree with the proposed $3\Theta$ correction.

For $d=3$, the exact computation is

$$
(1,1,-11/2,37/6,-71/24)e^{2\Theta}
=(1,3,-3/2,-3/2,3/8)=\alpha+3\beta.
$$

## Additional checks

The program calculates exact rational results for each transformation in this
list:

- theta twists
- dualization
- shifts and cone reversal
- translation
- degree-zero line bundles
- the Todd/Mukai convention
- multiplication pullback
- parameter relabeling
- component counts
- point-normalization counts

The degree-two equation gives two nondual line-twist candidates. Only the
candidate with the additional $+2\Theta$ twist also satisfies the degree-three
and degree-four equations. Dualization gives the related
$\alpha-3\beta$ solution.

The correction does not change this normalized characteristic class:

$$
\kappa(E)=\mathrm{ch}(E)
\exp\left(-c_1(E)/\mathrm{rank}(E)\right).
$$

The Hodge-side contraction has the same rank before and after the correction.
In both cases, the rank is $12$, and the kernel dimension is $16$. Only the
mixed witness changes. For $d=3$, the witness is

$$
C_{ij}-3A_{ij}.
$$

The corrected class produces this secant-plane form. This result is a check of
one consequence of the correction. It does not supply the Atiyah/evaluation
map that the construction still needs.

## The separate $\alpha+\beta$ construction

The repaired Section 12 theta-cycle has character $\alpha+3\beta$. It does not
have the $\alpha+\beta$ character of the older genus-four Jacobian construction.
That construction appears in Example 8.2.3 of
[arXiv:2502.03415v2](https://arxiv.org/abs/2502.03415).

No exact derived autoequivalence of $D^b(X)$ can identify these two classes.
Derived autoequivalences preserve the Euler pairing. Their self-Euler pairings
are

$$
\chi(\alpha+3\beta,\alpha+3\beta)=8d(d+9),\qquad
\chi(\alpha+\beta,\alpha+\beta)=8d(d+1).
$$

At $d=3$, these values are $288$ and $96$. The older $W_2$ construction
realizes $\alpha+\beta$ exactly. It remains a separate construction.

## Reproduce the results

These checks require only the Python standard library. Run these commands:

```bash
PYTHONPATH=src python3 -m thetatwist --d 3
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

Use this command to rebuild the committed JSON result and its SHA-256
manifest:

```bash
PYTHONPATH=src python3 -m thetatwist --d 3 --write-results results
```

[PROOF.md](PROOF.md) contains the full coefficient derivation.
[results/d3_certificate.json](results/d3_certificate.json) contains the
generated machine result.

The experimental continuation is in
[EVALUATION_CAMPAIGN.md](EVALUATION_CAMPAIGN.md). It uses only the corrected
Section 12 object. It gives the exact theta-incidence complex and reduces the
open Atiyah calculation from 16 kernel directions to six mixed directions.

## Citation and license

[CITATION.cff](CITATION.cff) contains the citation metadata. The repository
uses the [MIT License](LICENSE).

## Scope

This repository verifies a Chern-character identity and its immediate
Hodge-contraction consequences. It does not prove semiregularity. It does not
give the Atiyah/evaluation map that the construction still needs. It does not
prove a new case of the Hodge conjecture.

Cyclic divisor indexing does not define an effective cyclic group action
without separate geometric data.

Brian Ross prepared this work independently. He is a graduate student at the
University of Oklahoma. This repository is not an institutional publication.
It does not make a claim about the intent of the paper's author. It presents a
proposed correction. It proves uniqueness among integral theta twists within
the stated cyclic presentation.
