"""Symbolic model inversion with SymPy, called from R via reticulate."""

import sympy
from sympy import Symbol, sympify, solve, diff, simplify, oo


def invert_model(model_str, bio_par, all_symbols):
    """Solve ``model == value`` for the biological parameter.

    Parameters
    ----------
    model_str : str
        Model as given by the user, in SymPy syntax, e.g. ``"yi / sj"``.
    bio_par : str
        Name of the parameter to solve for.
    all_symbols : list[str]
        Names of all parameter symbols in ``model_str``.

    Returns
    -------
    dict or None
        ``{"inverse": str, "inverse_deriv": str}`` in R syntax, or ``None``
        if no verified closed-form inverse exists.
    """
    sym_dict = {}
    for name in list(all_symbols) + ["value"]:
        sym_dict[name] = Symbol(name, positive=True)

    # Same logarithms as stats::D() supports; unregistered names would parse
    # as undefined functions that solve() cannot invert.
    sym_dict["log2"] = lambda a: sympy.log(a, 2)
    sym_dict["log10"] = lambda a: sympy.log(a, 10)

    try:
        model_expr = sympify(model_str, locals=sym_dict)
    except Exception:
        return None

    equation = model_expr - sym_dict["value"]
    target = sym_dict[bio_par]

    try:
        solutions = solve(equation, target)
    except Exception:
        return None

    if not solutions:
        return None

    inverse_expr = _pick_positive_solution(solutions, sym_dict)
    if inverse_expr is None:
        return None

    # Inverses with transcendental functions of value are rejected.
    if _contains_transcendental(inverse_expr, sym_dict["value"]):
        return None

    try:
        inverse_deriv = diff(inverse_expr, sym_dict["value"])
        inverse_deriv = simplify(inverse_deriv)
    except Exception:
        return None

    try:
        check = simplify(equation.subs(target, inverse_expr))
        if check != 0:
            return None
    except Exception:
        return None

    inverse_str = _sympy_to_r(inverse_expr)
    deriv_str = _sympy_to_r(inverse_deriv)

    return {"inverse": inverse_str, "inverse_deriv": deriv_str}


def _pick_positive_solution(solutions, sym_dict):
    """Return the only solution that is positive for positive symbols, else None."""
    if len(solutions) == 1:
        return solutions[0]

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
    return None


def _contains_transcendental(expr, value_sym):
    """Whether expr applies a transcendental function to value_sym."""
    transcendental_names = {
        "sin", "cos", "tan", "asin", "acos", "atan", "atan2",
        "sinh", "cosh", "tanh", "asinh", "acosh", "atanh",
        "LambertW",
    }
    for func in expr.atoms(sympy.Function):
        if str(type(func).__name__) in transcendental_names:
            if value_sym in func.free_symbols:
                return True
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
