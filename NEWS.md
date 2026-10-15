# blotIt 0.4.1

* `alignReplicates()` gains `tGroup` for `likelihood = "t"`: the points of a
  group (e.g. a gel) follow a joint multivariate t distribution and share one
  scale weight, so a noisy group is down-weighted as a whole.

# blotIt 0.4.0

* `alignReplicates()` gains `likelihood = c("normal", "t")` and `df`. With
  `"t"`, each data point follows a Student-t distribution with scale
  `errorModel`; uncertainties use its expected Fisher information. The default
  `"normal"` gives the same results as before.
* The result has `likelihood` and `likelihoodDf`; `llrTest()` refuses to
  compare fits with different likelihoods.
* Passing a previous `alignReplicates()` result as `data` refits its data and
  starts at its estimates; `parinit` takes precedence.
* Output tables no longer contain the internal column `"1"`.
* `llrTest()` checked the nesting of error effects against the scaling
  effects.
* `plotIt()` no longer fails for `plotLine = "original"` with points from
  `"scaled"` or `"prediction"`, and uses current ggplot2 syntax
  (requires ggplot2 >= 3.4.0).

# blotIt 0.3.0

* `alignReplicates()` gains `parinit` (user start values per parameter or
  level), `multistartSd` (spread of the multistart perturbations) and `cores`
  (multistart fits in parallel, forked).
* `error` accepts a list of formulas, one per error parameter group.
* Parameters without influence on the likelihood are held at their start value
  and reported as `NA` instead of failing on a singular Hessian.
* The result has a `multistart` table; `plotWaterfall()` plots it.
* Parameter names containing digits (e.g. `sigma0`) are classified correctly.
* The default call returns the same results as 0.2.0.
