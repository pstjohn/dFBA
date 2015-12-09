cdef class ReactionBoundsFunction:
    cdef int evaluate(self, realtype* y, double [:] lb_e, double [:] ub_e):
        """ A user-defined function which updates the lower and upper bound of the
        metabolic exchange fluxes.

        The current external concentrations, y, are passed as inputs. lb_e and ub_e
        are passed with their current (default) values, so only active bounds need
        to be changed.
        """

        return 0
