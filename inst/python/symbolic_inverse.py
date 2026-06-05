"""
Symbolic model inversion using SymPy.

Called from R via reticulate.  The module lives in inst/python/ and is
loaded with import_from_path() + delay_load following the CppODE pattern.
"""

import sympy
from sympy import Symbol, sympify, solve, diff, simplify, oo


def invert_model(model_str, bio_par, all_symbols):
    """Solve  model == value  for bio_par and return inverse + derivative.

    Parameters
    ----------
    model_str : str
        The model equation as typed by the user, e.g. "yi / sj",
        "yi * sj + offset".  Must be valid SymPy syntax.
    bio_par : str
        Name of the biological parameter to solve for, e.g. "yi".
    all_symbols : list[str]
        Names of *all* parameter symbols in ``model_str``.

    Returns
    -------
    dict or None
        On success: ``{"inverse": <R string>, "inverse_deriv": <R string>}``.
        ``None`` when SymPy cannot solve analytically or verification fails.
    """
    # ---- declare symbols (all assumed positive & real) --------------------
    sym_dict = {}
    for name in list(all_symbols) + ["value"]:
        sym_dict[name] = Symbol(name, positive=True)

    # ---- register base-2 / base-10 logarithms -----------------------------
    # SymPy only knows the natural ``log``; ``log2``/``log10`` would otherwise
    # be parsed as undefined functions and ``solve`` could not invert them.
    # stats::D() in the forward model supports exactly these two (but not the
    # two-argument ``log(x, base)`` form), so we mirror that here.
    sym_dict["log2"] = lambda a: sympy.log(a, 2)
    sym_dict["log10"] = lambda a: sympy.log(a, 10)

    # ---- parse model & build equation: model - value = 0 -----------------
    try:
        model_expr = sympify(model_str, locals=sym_dict)
    except Exception:
        return None

    equation = model_expr - sym_dict["value"]
    target = sym_dict[bio_par]

    # ---- solve -----------------------------------------------------------
    try:
        solutions = solve(equation, target)
    except Exception:
        return None

    if not solutions:
        return None

    # ---- pick the right solution -----------------------------------------
    # For unique solutions, just take it.
    # For multiple solutions, prefer the one that is positive when all
    # symbols are positive (which is the physical expectation in blotIt).
    inverse_expr = _pick_positive_solution(solutions, sym_dict)
    if inverse_expr is None:
        return None

    # ---- reject transcendental / non-algebraic inverses ------------------
    # If the inverse contains trig, inverse-trig, or other transcendental
    # functions of 'value', it is not a clean algebraic inverse.
    if _contains_transcendental(inverse_expr, sym_dict["value"]):
        return None

    # ---- derivative d(inverse)/d(value) ----------------------------------
    try:
        inverse_deriv = diff(inverse_expr, sym_dict["value"])
        inverse_deriv = simplify(inverse_deriv)
    except Exception:
        return None

    # ---- verify: substitute back, should simplify to 0 -------------------
    try:
        check = simplify(equation.subs(target, inverse_expr))
        if check != 0:
            return None
    except Exception:
        return None

    # ---- convert to R syntax strings ------------------------------------
    inverse_str = _sympy_to_r(inverse_expr)
    deriv_str = _sympy_to_r(inverse_deriv)

    return {"inverse": inverse_str, "inverse_deriv": deriv_str}


def _pick_positive_solution(solutions, sym_dict):
    """From a list of SymPy solutions, pick the one that evaluates
    positive when all symbols are positive.  Returns None if ambiguous."""
    if len(solutions) == 1:
        return solutions[0]

    # Try substituting large positive values for all symbols
    test_subs = {s: 2 for s in sym_dict.values()}
    positive_sols = []
    for sol in solutions:
        try:
            val = float(sol.subs(test_subs))
            if val > 0:
                positive_sols.append(sol)
        except Exception:
            continue

    if len(positive_sols) == 1:
        return positive_sols[0]

    # Ambiguous or no positive solution
    return None


def _contains_transcendental(expr, value_sym):
    """Check whether expr contains transcendental functions of value_sym."""
    # Walk the expression tree for functions applied to terms involving value
    transcendental_names = {
        "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
        "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
        "LambertW",
    }
    for func in expr.atoms(sympy.Function):
        if str(type(func).__name__) in transcendental_names:
            if value_sym in func.free_symbols:
                return True
    # Also check for applied unevaluated functions
    for sub in sympy.preorder_traversal(expr):
        if isinstance(sub, sympy.Function):
            fname = type(sub).__name__
            if fname in transcendental_names and value_sym in sub.free_symbols:
                return True
    return False


def _sympy_to_r(expr):
    """Convert a SymPy expression to an R-compatible string."""
    s = str(expr)
    s = s.replace("**", "^")
    return s
