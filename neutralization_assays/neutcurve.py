"""
=============
neutcurve
=============
Module for fitting and analyzing neutralization curves.
"""


def fourParamLogistic(c, ic50, slope, top, bottom):
    """4-parameter logistic neutralization curve.

    Evaluates :math:`f(c) = b + \\frac{t - b}{1 + (c/i)^s}`

    In the argument definitions below, we interpret the
    meaning of this curve in the case where :math:`s > 0`,
    in which case :math:`f(c)` can be thought of as
    the fraction surviving at concentration `c`.

    Args:
        `c` (float)
            Concentration of antibody
        `ic50` (float)
            Inhibitory concentration 50%, :math:`i` in equation.
        `slope` (float)
            Hill slope of curve, :math:`s` in equation.
        `bottom` (float)
            Bottom of curve (value as :math:`c` gets large),
            :math:`b` in equation.
        `top` (float)
            Top of curve (value as :math:`c` get small),
            :math:`t` in equation.

    Returns:
        A float giving :math:`f(c)`, which is interpreted
        as the fraction surviving at concentration `c`.
    """
