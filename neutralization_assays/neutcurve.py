"""
=============
neutcurve
=============
Module for fitting and analyzing neutralization curves.
"""


import functools
import scipy
import scipy.optimize


class fourParamLogistic:
    """4-parameter logistic neutralization curve.

    Evaluates :math:`f(c) = b + \\frac{t - b}{1 + (c/i)^s}`

    When :math:`s > 0`, then :math:`f(c)`, then :math:`f(c)`
    is the fraction surviving at concentration :math:`c`.

    Each object has the following attributed:
        `ic50` (float)
            Inhibitory concentration 50%, :math:`i` in equation
            above
        `slope` (float)
            Hill slope of curve, :math:`s` in equation above.
        `bottom` (float)
            Bottom of curve (value as :math:`c` gets large),
            :math:`b` in equation above.
        `top` (float)
            Top of curve (value as :math:`c` get small),
            :math:`t` in equation above.

    See `__init__` method for how to initialize an object.

    To use an object to evaluate neutralization values:

    >>> top = 1
    >>> bottom = 0
    >>> slope = 5
    >>> ic50 = 0.03
    >>> neut = fourParamLogistic(ic50, slope, top, bottom)
    >>> scipy.allclose(0.5, neut.fracsurvive(ic50))
    True
    >>> scipy.allclose(bottom, neut.fracsurvive(1e3))
    True
    >>> scipy.allclose(top, neut.fracsurvive(1e-3))
    True
    """


    def __init__(self, ic50=0, slope=1, top=1, bottom=0):
        """Initialize a new object."""
        self.ic50 = ic50
        self.slope = slope
        self.top = top
        self.bottom = bottom


    def fracsurvive(self, c):
        """Evaluates fraction surviving at concentration `c`."""
        return self._fracsurvive_allargs(c, self.ic50,
                self.slope, self.top, self.bottom)


    def _fracsurvive_allargs(self, c, ic50, slope, top, bottom):
        """Evaluates fraction surviving with called params.

        Differs from `fracsurvive` in that it does **not** use
        the current object values for parameters.
        """
        return bottom + (top - bottom) / (1 + (c / ic50)**slope)


    def fitcurve(self, cs, fs, fixtop=False, fixbottom=False):
        """Fits neutralization curve parameters to data.

        Using the provided data, fits values for the 
        neutralization curve attributes.

        Args:
            `cs` (array-like)
                List of concentrations for which we have
                measurements to fit.
            `fs` (array-like)
                List of same length as `cs`, with `fs[i]`
                giving fraction surviving at `cs[i]`.
            `fixtop` (bool)
                If `True`, keep the current value of the
                `top` attribute and do not fit.
            `fixbottom` (bool)
                If `True`, keep the current value of the
                `bottom` attribute and do not fit.
        """
        # get data into arrays sorted by concentration
        cs = scipy.array(cs)
        fs = scipy.array(fs)
        fs = fs[cs.argsort()]
        cs = cs[cs.argsort()]

        # make initial guess for slope have the right sign
        if fs[0] > fs[-1]:
            slopeguess = 1.0
        else:
            slopeguess = -1.0

        # guess ic50
        midval = (fs.max() - f.min()) / 2.0

        raise RuntimeError("need to deal with midpoint versus IC50")

        if fixtop and fixbottom:
            func = functools.partial(self._fracsurvive_allargs,
                    top=self.top, bottom=self.bottom)
        elif fixtop:
            func = functools.partial(self._fracsurvive_allargs,
                    top=self.top)
        elif fixbottom:
            func = functools.partial(self._fracsurvive_allargs,
                    bottom=self.bottom)
        else:
            func = self._fracsurvive_allargs

        # some intelligent guesses for initial values

        scipy.optimize.curve_fit(
                func,
                cs,
                fs,
                initguess
                )




if __name__ == '__main__':
    import doctest
    doctest.testmod()
