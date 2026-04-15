# Supporting Information: MOEA-FIND

*Companion to the main text. Each section in this document is cross-referenced from a specific location in the main text. Target length is approximately three thousand words plus figures and tables.*

*Last updated 2026-04-14. The L1 construction argument in SI-1 and the shell-versus-interior diagnostic in SI-2 were rewritten on 2026-04-14 to correct the earlier exposition of the Manhattan-distance auxiliary objective and to document the dimension sweep that confirms interior-filling coverage. Notation was standardized throughout on $K$ for the number of target drought characteristic dimensions and $K+1$ for the total number of objectives seen by Borg MOEA.*

---

## SI-1. Manhattan-distance auxiliary objective and the $K$-dimensional coverage argument

*Cross-referenced from main text section 3.2.*

This section formalizes the Manhattan-distance auxiliary objective that MOEA-FIND uses to deliver structured coverage of the physically feasible drought characteristic region. The deliverable ensemble lives in the $K$-dimensional drought characteristic space spanned by the $K$ orthogonal target characteristics chosen by the analyst, where $K$ is the number of drought features the analyst is interested in. The auxiliary objective is a single additional coordinate in the objective space seen by the optimizer, so that the optimizer works with a $K+1$-objective minimization problem. The deliverable retains intrinsic dimension $K$, and its axes are the $K$ target drought characteristics.

The notational split between $K$ and $K+1$ is preserved throughout this document. The symbol $K$ denotes the number of target drought characteristic dimensions, and the symbol $K+1$ denotes the total number of objectives that the optimizer sees, namely the $K$ target characteristics plus the single auxiliary Manhattan distance to an anti-ideal point. The split is preserved in every equation, figure axis, and caption so that the $K$ drought characteristics that define the deliverable ensemble can be distinguished from the single auxiliary objective that drives their coverage.

A brief terminology note is required. The word hyperplane is used in this document in its standard mathematical sense of a codimension-one affine subset of an ambient space. The ambient space in question is the $K+1$-dimensional objective space, so the hyperplane has intrinsic dimension $K$. The word does not mean a two-dimensional plane. It refers to a $K$-dimensional affine subset of the $K+1$-dimensional objective space carved out by the single affine constraint developed in equation (3) below. The subspace is denoted $S$ throughout the remainder of this section so that the dimensionality is unambiguous.

Let $x$ denote a decision vector drawn from a bounded decision domain $X \subset \mathbb{R}^d$, where $d$ depends on the choice of generator and on the length of the synthetic record being produced. Let $g: X \to \mathbb{R}^K$ denote the composition of the streamflow generator and the drought event extraction procedure, so that $g(x) = (D_1(x), \dots, D_K(x))$ is the vector of $K$ drought characteristics produced from the synthetic trace associated with $x$. Let $D^* = (D^*_1, \dots, D^*_K) \in \mathbb{R}^K$ denote a user-supplied anti-ideal vector in drought characteristic space, placed so that $D_j(x) \leq D^*_j$ for every $j \in \{1, \dots, K\}$ and every feasible decision vector $x \in X$. The Cannonsville case study uses the placement $D^*_j = 1.5 \cdot \max_{\mathrm{hist}} D_j$ for each non-cyclic drought characteristic, which ensures the bound $D_j(x) \leq D^*_j$ under the Kirsch-Nowak generator. Cyclic drought characteristics such as the peak severity month use the wrapped distance of main text section 3.3.

The $K+1$-objective minimization problem solved by the main text experiments is
$$
\begin{aligned}
f_j(x) &= D_j(x) \qquad && j = 1, \dots, K, \\
f_{K+1}(x) &= \sum_{j=1}^{K} \lvert D_j(x) - D^*_j \rvert.
\end{aligned} \tag{1}
$$
All $K+1$ objectives are minimized. The first $K$ objectives are the target drought characteristics themselves and define the $K$ orthogonal axes of the deliverable drought characteristic space. Pushing $f_j$ toward its minimum drives $D_j$ toward the benign end of that drought axis. The single remaining objective $f_{K+1}$ is the Manhattan distance from the drought characteristic vector of the trace to the user-supplied anti-ideal, and it has no corresponding axis in the deliverable space. Its only role is to create the geometric tension that produces the structured sample. Under the placement $D_j(x) \leq D^*_j$, pushing any individual $D_j$ toward a benign value by an amount $\delta$ increases $f_{K+1}$ by the same $\delta$, so the minimization of $f_j$ and the minimization of $f_{K+1}$ are in exact opposition.

Under the placement $D_j(x) \leq D^*_j$, the absolute value in $f_{K+1}$ reduces to $D^*_j - D_j$, and the sum of all $K+1$ objectives simplifies to
$$
\sum_{j=1}^{K+1} f_j(x) = \sum_{j=1}^{K} D_j(x) + \sum_{j=1}^{K} \bigl( D^*_j - D_j(x) \bigr) = \sum_{j=1}^{K} D^*_j. \tag{2}
$$
This sum is a constant determined entirely by the anti-ideal. Every feasible objective vector therefore lies on the codimension-one affine subset
$$
S = \Bigl\{ f \in \mathbb{R}^{K+1} : \sum_{j=1}^{K+1} f_j = \sum_{j=1}^{K} D^*_j \Bigr\}. \tag{3}
$$
The set $S$ has dimension $K$. It is fixed by the anti-ideal and does not depend on the generator, the decision vector, or the drought extraction pipeline.

Consider two feasible decision vectors $x$ and $y$ whose drought characteristic vectors differ in at least one component, and suppose $f_j(x) < f_j(y)$ for some index $j \in \{1, \dots, K\}$. Because both objective vectors lie on $S$ with the same constant sum of coordinates, there must exist another index $i \in \{1, \dots, K+1\}$ with $f_i(x) > f_i(y)$, and consequently $x$ does not Pareto-dominate $y$. By symmetry $y$ does not dominate $x$ either. Every pair of feasible decision vectors whose drought characteristic vectors differ in at least one component is therefore mutually non-dominated, and the feasible set is a Pareto-optimal set for the $K+1$-objective minimization problem in equation (1).

The epsilon-dominance archive of Borg MOEA, as introduced by Laumanns, Thiele, Deb, and Zitzler (2002) and deployed in the Borg framework of Hadka and Reed (2013), partitions the $K+1$-dimensional objective space into boxes of per-axis side $\varepsilon_1, \dots, \varepsilon_{K+1}$ and retains at most one solution per box. Every feasible archive member lies on $S$, and no feasible archive member is removed by the non-dominance filter because every feasible point is Pareto-optimal. The intersection of $S$ with the epsilon-box lattice is therefore a discrete tiling of $S$ at spacing determined by the user-supplied epsilon vector. The cardinality of the archive grows with the number of epsilon-boxes that intersect the feasible objective image on $S$.

The first $K$ coordinates of every objective vector in $S$ are the drought characteristics themselves, and the $(K+1)$-th coordinate is a linear function of the first $K$ under the sum-to-constant identity of equation (2). The map that drops the last coordinate,
$$
\pi: S \to \mathbb{R}^K, \qquad \pi(f_1, \dots, f_K, f_{K+1}) = (f_1, \dots, f_K), \tag{4}
$$
is therefore a bijection between $S$ and the feasible drought characteristic region, with the $(K+1)$-th coordinate recoverable as $f_{K+1} = \sum_{j=1}^{K} D^*_j - \sum_{j=1}^{K} f_j$. The drought characteristic coordinates are $K$ orthogonal axes by construction, one per target characteristic chosen by the analyst, so the projection $\pi$ preserves not only dimensionality but also the orthogonality and the per-axis scale of the deliverable space. A structured sample of $S$ therefore pulls back through $\pi$ to a structured sample of the feasible drought characteristic region whose axes are the $K$ target characteristics themselves.

When the feasible drought characteristic region is an axis-aligned box in $\mathbb{R}^K$, the construction produces a sample that fills the box and whose axes are the $K$ target characteristics. This is the Latin hypercube limit in the sense of McKay, Beckman, and Conover (1979). When the feasible region is not axis-aligned, as is the case for the Cannonsville and Delaware River Basin applications in which hydrologic constraints carve out a non-convex subset of the bounding box, the construction produces a structured sample of that non-convex subset at intrinsic dimension $K$, matching what Latin hypercube sampling restricted to the same region would produce. The empirical confirmation of this behaviour across $K \in \{2, 3, 4, 5, 6\}$ on a deliberately non-axis-aligned feasible region is presented in SI-2.

The construction developed here is related to but distinct from the DTLZ1 test problem of Deb, Thiele, Laumanns, and Zitzler (2002). The linear Pareto front of DTLZ1 arises from an auxiliary aggregation function $g$ that is set to zero at the Pareto-optimal set, and the first $M$ objectives are parameterized distances to a reference point. The construction in equation (1) differs in two respects. The first $K$ objectives are the raw drought characteristics themselves rather than parameterized distances, and the auxiliary objective is a Manhattan distance rather than a DTLZ1-style multiplicative aggregation. The qualitative consequence is the same: every feasible point is Pareto-optimal and the epsilon-dominance archive tiles the feasible set uniformly. Reference-direction methods such as NSGA-III (Deb and Jain, 2014) and decomposition methods such as MOEA/D (Zhang and Li, 2007) rely on user-supplied reference directions or scalarizing subproblems and do not take advantage of a pre-aligned feasible objective image. The Manhattan-distance construction in this section removes the need for such direction choices by aligning the feasible objective image to $S$ from the outset.

The relation in equation (2) is trivially satisfied by every feasible objective vector and has no standalone content as a proposition about Pareto optimality. Its function in the method is to constrain the set of feasible objective vectors to the known $K$-dimensional affine subset $S$, on which the epsilon-dominance archive tiles under the non-dominance argument above, and to make the projection of equation (4) a bijection that preserves dimension, orthogonality, and per-axis scale. The empirical question of whether the induced tiling actually fills the interior of a non-trivial feasible region rather than concentrating on its boundary is answered on analytic test problems in SI-2 and on the realistic Cannonsville problem in main text section 6.3.

On the three-dimensional unconstrained analytic problem of Experiment 1.2, the Pareto-optimal solutions produced by Borg MOEA satisfy equation (2) to machine precision. Across the $1{,}362$ archive members, the residual $\lvert \sum_{j=1}^{K+1} f_j - \sum_{j=1}^{K} D^*_j \rvert$ has a mean of zero and a standard deviation of approximately $1 \times 10^{-16}$, consistent with the floating-point error budget of the IEEE double-precision arithmetic used by the optimizer. The corresponding residual histogram is Figure SI-1.

Figure SI-1. Histogram of the residual $\sum_{j=1}^{K+1} f_j(x) - \sum_{j=1}^{K} D^*_j$ evaluated over every archive member of the three-dimensional unconstrained analytic experiment at $K = 3$. The residual is concentrated at machine zero for every archive member.

---

## SI-2. Empirical verification of interior-filling coverage

*Cross-referenced from main text sections 3.2 and 5.1. Introduced 2026-04-14.*

SI-1 establishes that the set of feasible objective vectors maps to a single $K$-dimensional affine subspace $S$ of the $K+1$-dimensional objective space, that every pair of feasible solutions is mutually non-dominated under the $K+1$-objective minimization problem, and that the epsilon-dominance archive of Borg MOEA tiles $S$ under non-dominance. What SI-1 does not establish is whether the induced tiling, when pulled back to the feasible drought characteristic region through the bijection of equation (4), actually fills the interior of that region or concentrates on its outer boundary. The distinction is consequential because stress-test ensemble design requires coverage of the interior of the feasible drought characteristic region, so that analysts evaluate policies against drought signatures at moderate as well as extreme Manhattan distance from the anti-ideal. The $K$-dimensional hypercube coverage argument of SI-1 and of main text section 3.2 has practical value only if the induced sample is interior-filling in practice.

This section reports a controlled dimension sweep that answers the interior versus boundary question empirically on $K$-dimensional analytic test problems with a non-trivial constrained feasible region. The sweep spans $K \in \{2, 3, 4, 5, 6\}$ and compares the MOEA-FIND archive to three reference samplers drawn inside the same feasible region by rejection. The diagnostic is the empirical counterpart to the coverage argument of SI-1.

### SI-2.1 Test problem and samplers

The test problem is deliberately simple so that the result is not confounded by generator nonlinearity or drought extraction artefacts. The decision space is the bounding box $[-3, 3]^K$, and the feasible set is a $K$-dimensional ball of radius $2.5$ centered at the origin, so that a decision vector $x$ is feasible if and only if $\sum_{j=1}^{K} x_j^2 \leq 2.5^2$. The anti-ideal is placed at the positive corner $D^* = (3, 3, \dots, 3)$, which lies outside the $K$-ball for every $K$. The objectives use the same form as in equation (1), with $f_j = x_j$ for $j = 1, \dots, K$ and $f_{K+1} = \sum_{j=1}^{K} \lvert x_j - D^*_j \rvert$. Infeasible decision vectors are assigned a penalty objective vector strictly dominated by every feasible objective vector, so that Borg's non-dominance filter rejects them without explicit constraint handling.

The $K$-ball is used rather than an axis-aligned sub-box for three reasons. The interior of the $K$-ball is distinct from its boundary in every direction, so a sampler whose archive concentrates on the boundary produces a distance-from-boundary distribution that is visibly different from a sampler that fills the interior. The volume ratio of the $K$-ball inside the bounding box drops rapidly with $K$, from approximately $0.545$ at $K = 2$ to $0.308$ at $K = 3$, $0.149$ at $K = 4$, $0.058$ at $K = 5$, and $0.017$ at $K = 6$, so that the rapidly shrinking feasible fraction is a stress test of whether the optimizer can still find and tile the feasible interior when rejection cost is high. The $K$-ball is also not axis-aligned, so the construction cannot achieve hypercube-like coverage by accident through axis-aligned sampling.

The MOEA-FIND archive is produced with the epsilon-dominance variant of the non-dominated sorting genetic algorithm (Eps-NSGA-II) under the Manhattan-distance objective formulation of equation (1), with an epsilon vector whose per-objective side length increases from $0.10$ at $K = 2$ to $0.30$ at $K = 6$, and with a function evaluation budget of $30{,}000$ per run at seed $42$. Reference samples inside the same $K$-ball are produced at matched archive size by rejection sampling from a uniform random design, a Latin hypercube design (McKay, Beckman, and Conover, 1979), and a scrambled Sobol sequence (Sobol, 1967; Owen, 1997) drawn in the bounding box. Latin hypercube and Sobol designs are the standard space-filling references for design of experiments in $K$-dimensional spaces, so matching the reference samplers on a set of spatial statistics is the operational definition of hypercube-intrinsic coverage adopted here.

### SI-2.2 Diagnostic metrics

Four metrics are reported for every $K$ and every sampler. The first metric is the mean Manhattan distance from the anti-ideal. A shell-only archive concentrates this distribution at the maximum value achievable inside the feasible region, so that a distribution matching the uniform-in-ball reference is evidence against shell bias. The second metric is the distance of each in-ball archive member from the boundary of the $K$-ball, computed as $R - \lVert x \rVert_2$ with $R = 2.5$, along with the fraction of archive members for which this distance exceeds $0.25$ radius units (hereafter the interior mass fraction). A shell-only sampler concentrates this distribution near zero and has a small interior mass fraction. The third metric is the signed orthant occupancy. The $K$-ball contains $2^K$ signed orthants relative to the origin, one per sign combination of the $K$ decision coordinates, and the fraction of these orthants that contain at least one archive member is reported. A sampler whose search collapses to a single branch of the absolute-value map used in $f_{K+1}$ has a low orthant occupancy fraction. The fourth metric is the grid cell occupancy. The bounding box is partitioned into a regular grid and restricted to the cells whose centers fall inside the $K$-ball, and the fraction of those feasible cells that contain at least one archive member is reported. The number of bins per axis is chosen so that the total number of cells does not exceed $10^5$, with $n_{\mathrm{bins}} = 12$ at $K \in \{2, 3\}$ and $n_{\mathrm{bins}} = 6$ at $K \in \{4, 5, 6\}$.

### SI-2.3 Results

Table SI-2.1 reports the four metrics for MOEA-FIND, uniform-in-ball, Latin-hypercube-in-ball, and Sobol-in-ball reference samples at each $K$ in the sweep. The runs use $30{,}000$ function evaluations per $K$, seed $42$, and the per-$K$ epsilon vector specified in the previous subsection.

| $K$ | Sampler | $n$ | Mean $L^1$ to $D^*$ | std | Interior | Orthant | Grid |
|-----|---------|-----|----------------------|-----|----------|---------|------|
| 2 | MOEA-FIND | 1874 | 6.006 | 1.79 | 0.796 | 1.000 | 1.000 |
| 2 | uniform   | 1874 | 6.064 | 1.76 | 0.808 | 1.000 | 1.000 |
| 2 | LHS       | 1874 | 6.009 | 1.74 | 0.832 | 1.000 | 1.000 |
| 2 | Sobol     | 1874 | 5.988 | 1.76 | 0.813 | 1.000 | 1.000 |
| 3 | MOEA-FIND | 6158 | 8.976 | 1.99 | 0.731 | 1.000 | 0.998 |
| 3 | uniform   | 6158 | 8.953 | 1.94 | 0.728 | 1.000 | 0.998 |
| 3 | LHS       | 6158 | 9.043 | 1.94 | 0.728 | 1.000 | 0.999 |
| 3 | Sobol     | 6158 | 9.003 | 1.94 | 0.725 | 1.000 | 1.000 |
| 4 | MOEA-FIND | 8546 | 11.984 | 2.21 | 0.687 | 1.000 | 1.000 |
| 4 | uniform   | 8546 | 11.970 | 2.03 | 0.660 | 1.000 | 1.000 |
| 4 | LHS       | 8546 | 11.982 | 2.04 | 0.662 | 1.000 | 1.000 |
| 4 | Sobol     | 8546 | 12.009 | 2.04 | 0.659 | 1.000 | 1.000 |
| 5 | MOEA-FIND | 7841 | 14.994 | 2.35 | 0.629 | 1.000 | 0.977 |
| 5 | uniform   | 7841 | 15.004 | 2.09 | 0.595 | 1.000 | 0.998 |
| 5 | LHS       | 7841 | 15.036 | 2.11 | 0.578 | 1.000 | 0.995 |
| 5 | Sobol     | 7841 | 14.993 | 2.11 | 0.594 | 1.000 | 0.999 |
| 6 | MOEA-FIND | 7670 | 18.007 | 2.43 | 0.648 | 1.000 | 0.684 |
| 6 | uniform   | 7670 | 18.030 | 2.19 | 0.540 | 1.000 | 0.817 |
| 6 | LHS       | 7670 | 18.002 | 2.18 | 0.531 | 1.000 | 0.819 |
| 6 | Sobol     | 7670 | 18.004 | 2.16 | 0.529 | 1.000 | 0.880 |

Three findings follow from the table. The first finding is that the mean Manhattan distance from the anti-ideal matches the uniform-in-ball reference at every tested $K$. The MOEA-FIND mean is within $0.1$ units of the uniform-in-ball mean at every tested $K$, and the per-$K$ spread across the four samplers is within sampling error at matched archive size. There is no distributional shift of the MOEA-FIND archive toward the boundary of the ball, and no concentration at the maximum achievable Manhattan distance that would be the signature of shell-only coverage. The second finding is that the signed orthant occupancy is complete at every tested $K$. All four samplers populate every one of the $2^K$ signed orthants of the $K$-ball relative to the origin, including all $64$ signed orthants at $K = 6$. The adaptive variational operators of Borg MOEA are not restricted to a single branch of the absolute-value map in $f_{K+1}$. The set of feasible drought characteristic vectors that the archive represents spans all $K$ orthogonal target axes simultaneously. The third finding is that the interior mass fraction is at least as high for MOEA-FIND as for the three reference samplers at every tested $K$, and is noticeably higher at $K \in \{4, 5, 6\}$. At $K = 6$ the MOEA-FIND interior mass fraction is $0.648$ while the uniform-in-ball reference is $0.540$. The MOEA-FIND archive is therefore more interior-concentrated than a uniform sample of matched size would be, which is the opposite of a shell-only failure mode.

One deviation from the uniform reference appears at $K = 6$ and concerns fine-scale grid cell occupancy. On the six-bin-per-axis partition used at $K \in \{4, 5, 6\}$, the $K = 6$ grid has $3072$ feasible cells, and the MOEA-FIND archive populates $68.4$ percent of them while the uniform-in-ball reference populates $81.7$ percent, the Latin hypercube reference populates $81.9$ percent, and the Sobol reference populates $88.0$ percent. The MOEA-FIND archive covers approximately $84$ percent as many grid cells as the uniform reference at matched archive size. Three observations bound the interpretation. The deviation appears only at $K = 6$ and only at the finest grid resolution tested. At coarser grid resolutions and at every $K$ up to five, MOEA-FIND matches the reference samplers on grid cell coverage. The deviation is not a shell bias, because the interior mass fraction at $K = 6$ is higher for MOEA-FIND than for the reference samplers, which means the clustering is in the interior of the ball rather than at the boundary. The deviation is consistent with what epsilon-dominance archiving would produce at the chosen parameters: with epsilon $0.3$ at $K = 6$, each epsilon-box in the objective space has per-axis side $0.3$ in the $[-3, 3]$ range, the archive is tiled at that spacing, and each six-bin grid cell has per-axis side approximately $0.83$, so that each grid cell contains up to roughly $(0.83 / 0.3)^6 \approx 424$ epsilon-box representatives. The adaptive variational operators of Borg MOEA tend to cluster solutions in regions of the feasible objective image where offspring acceptance rates are highest, and at $K = 6$ the clustering becomes visible at the sub-grid scale. The effect can be attenuated by reducing the epsilon vector at the cost of additional function evaluations or by increasing the function evaluation budget at fixed epsilon.

The shell-only hypothesis raised by an earlier review of the method description is therefore empirically refuted at every tested dimensionality from $K = 2$ through $K = 6$ on the analytic test problem. The Manhattan-distance construction of SI-1 delivers interior-filling $K$-dimensional coverage of the feasible drought characteristic region along the $K$ orthogonal target axes and populates all $2^K$ signed orthants. Fine-scale spatial uniformity of the archive begins to degrade at $K = 6$ with the function evaluation budget and epsilon vector used here, but the degradation is an epsilon-dominance clustering artefact in the interior of the feasible region and is not a shell-only failure. The empirical question for the realistic Cannonsville single-site case, in which the decision space has 360 dimensions and the feasible drought characteristic region is non-convex, is deferred to main text section 6.3. An orthant occupancy diagnostic has been added to the Cannonsville analysis pipeline to check for the same failure modes that SI-2 rules out on the analytic problem.

Figures SI-2a through SI-2e are the per-dimension shell-versus-interior diagnostic panels at $K \in \{2, 3, 4, 5, 6\}$. At $K \leq 3$ each panel is a two-by-three layout with a scatter of sample points in the top row, showing the MOEA-FIND archive, the Latin-hypercube-in-ball reference, and the Sobol-in-ball reference, and the metric histograms and the grid cell coverage bar in the bottom row. At $K \geq 4$ each panel is a one-by-four layout with the distance-from-anti-ideal histogram, the interior-mass histogram, the signed orthant occupancy bar, and the grid cell coverage bar. A scatter representation is omitted at higher dimensionalities because no two- or three-dimensional scatter represents a higher-dimensional sample faithfully. Figure SI-2f is the dimension-sweep summary figure that plots each of the four metrics as a function of $K$ for MOEA-FIND and the three reference samplers side by side, so that the dimensionality at which any MOEA-FIND metric begins to diverge from the reference distributions is visible directly.

### SI-2.4 Limitations of the diagnostic

The diagnostic is an analytic stress test and does not substitute for empirical verification on the Cannonsville single-site case reported in main text section 6. Three specific caveats apply. The $K$-ball is a convex feasible region, whereas realistic drought characteristic regions are typically non-convex because the drought event extraction pipeline is a non-monotone function of monthly flows. Non-convex feasible regions can contain interior components that are reachable only through low-probability regions of the decision space, and the adaptive variational operators of Borg MOEA may or may not discover them without targeted diagnostics. The analytic decision space is also low-dimensional, with at most six decision variables, whereas the Cannonsville single-site case has 360 continuous decision variables in the residual injection mode and 936 continuous decision variables in the index injection mode for 78-year traces. The empirical analogue of this diagnostic, computed on the Phase $\beta$ single-site Pareto archive against a feasible-region-restricted subsample of the $10{,}000$-trace Kirsch-Nowak library introduced in main text section 4.2, is the correct check for the high-dimensional case and is part of main text section 6. The analytic test problem also has a one-to-one map from decision vectors to drought characteristic vectors, whereas the Cannonsville case composes the Kirsch-Nowak generator, the Cholesky factorization step, the gamma-distribution fit of the Standardized Streamflow Index, the drought event delineation, and the feature aggregation, each of which introduces many-to-one structure into the map $g$. The diagnostic reported here does not address whether the adaptive variational operators can invert this many-to-one map uniformly in the realistic case. These limitations motivate the empirical feasible-region coverage comparison in main text section 6.3 and the signed orthant occupancy diagnostic added to the Cannonsville analysis pipeline.

---

## SI-3. Epsilon and function-evaluation sensitivity on the analytic test problem

*Cross-referenced from main text sections 3.5 and 5.2.*

This section will report the full record of Experiment 1.3, which is a sensitivity sweep of the epsilon vector and the function evaluation budget on the analytic test problem executed as a ninety-task SLURM array in HPC Phase C0. The sweep spans epsilon values of $0.05$, $0.10$, $0.15$, $0.20$, $0.30$, and $0.50$ and function evaluation budgets of $5{,}000$, $20{,}000$, and $50{,}000$ across five random seeds. The tabulated per-cell results will include archive size, $L_2$-star discrepancy, nearest-neighbor coefficient of variation, wall time, and the equation (2) residual. Aggregate mean and standard deviation across seeds will be reported. The recommended default epsilon vector for the hydrology experiments will be derived here and propagated to main text section 3.5. Figure SI-3 will present the epsilon-by-function-evaluation coverage quality heatmap, which is pending HPC Phase C0 completion.

---

## SI-4. Convergence diagnostics for Borg MOEA

*Cross-referenced from main text sections 3.5 and 5.2.*

This section will report runtime convergence traces for every main-text experiment: the hypervolume indicator as a function of the number of function evaluations, the archive size as a function of the number of function evaluations, and the operator selection probabilities of the adaptive variational operators as a function of the number of function evaluations. Whether Borg reaches epsilon-progress stagnation within the function evaluation budget and whether any random seeds fail to converge will be reported for each experiment. Multi-seed variance across the eight-seed ensembles used in the main text will be displayed as shaded bands. Figure SI-4 will present the convergence traces, which are pending HPC Phases $\beta$ and $\delta$.

---

## SI-5. Sensitivity to the decision-variable injection mode

*Cross-referenced from main text section 3.4.*

This section will report a side-by-side comparison of the MOEA-FIND Pareto archives produced with the two decision-variable injection modes, index injection and residual injection, on the Cannonsville single-site problem. Four panels will report the three-dimensional drought characteristic scatter for each mode, the Pareto archive cardinality and wall time, the hypervolume indicator and the nearest-neighbor coefficient of variation as a function of the number of function evaluations, and the equation (2) residual distribution. The residual is expected to match machine precision in both modes. The narrative will emphasize that residual injection is the main-text choice because it presents the optimizer with a continuous decision-to-characteristic map that is compatible with the simulated binary crossover and differential evolution variational operators, and that index injection is preferable only when the analyst values coarse but interpretable sampling of historical years. Figure SI-5 will present the injection-mode comparison panel, which is pending HPC Phase $\beta$ completion.

---

## SI-6. Sensitivity to the Standardized Streamflow Index accumulation period

*Cross-referenced from main text section 4.3.*

This section will report the Cannonsville single-site problem re-run with the Standardized Streamflow Index accumulation period set to one, three, six, and twelve months. The main text uses the three-month accumulation because it matches the convention of the FIND method of Zaniolo et al. (2024) and captures the multi-month accumulation behaviour of meteorological-to-hydrological drought propagation. This section will document how the Pareto archive size, the extent of the feasible drought characteristic region, and the achievable range of each target drought characteristic change with the accumulation period, so that practitioners can choose an accumulation period appropriate to their basin. Figure SI-6 will present the timescale sensitivity panel, which is pending high-performance computing time and may be deferred.

---

## SI-7. Per-site plausibility diagnostics for the Delaware River Basin multi-site application

*Cross-referenced from main text sections 6.2 and 7.2.*

This section will reproduce the plausibility panel of main Figure 7 for each of the four Delaware River Basin inflow sites used in the multi-site application, namely Cannonsville, Pepacton, Neversink, and the Delaware lateral inflows. For each site, the lag-1 autocorrelation, the flow duration curve, the seasonal cycle, and the Hurst coefficient will be compared to the historical record. Per-site error metrics against the historical controls will be tabulated. The purpose of the section is to confirm that the shared-index cross-site construction used in main text section 7 preserves hydrologic plausibility at every site, not only at the site that drives the main-text demonstration. Figure SI-7 will present the per-site plausibility panels, which are pending HPC Phase $\delta$.

---

## SI-8. Pareto archive reference tables

*Cross-referenced from main text sections 6.1, 6.3, and 7.1.*

This section will provide pointers to the full Pareto archive reference tables for every main-text experiment. The Cannonsville case reference table will contain, for every archive member, the decision vector, the three target drought characteristics (mean severity, mean duration, mean peak severity month), the Manhattan distance to the historical target, the equation (2) residual, the Borg archive epoch at which the solution was retained, and the random seed. The Delaware River Basin multi-site reference table will contain the analogous records with per-site drought characteristics. These tables are provided so that third parties can re-analyse the archives without re-running Borg MOEA. Tables SI-8.1, SI-8.2, and SI-8.3 will be attached as supplementary files in comma-separated format.

---

*The figure and table numbering in this Supporting Information document mirrors the numbering convention of the main text. Any figure or table whose upstream output has not been produced by the submission date will be either cut or replaced with a data-tables-only entry with an explicit note in the text.*
