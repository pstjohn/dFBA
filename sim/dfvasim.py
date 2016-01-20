import numpy as np
from .dfvasim_base import DFVA_Simulator_base


class DFVA_Simulator(object):

    def __init__(self, cobra_model, reaction_indicies, y0, bounds_func,
                 fva_tol=0.99):
        """ A class to handle the solution of the dynamic flux variablity
        analysis given in the passed model.
        
        cobra_model: a cobra.Model object
            A model containing the mass balance reactions and inner objective
            function (typically biomass formation)

        reaction_indicies: a numpy array
            A list of indicies corresponding to the external metabolites to be
            tracked

        y0: a numpy array
            An array of initial conditions for the specified external metabolites

        bounds_func: a ReactionBoundsFunction
            A function that specifies the uptake kinetics of the external
            metabolites. Specifically, this function adjusts the upper and
            lower bounds of the exchange rates as a function of current
            concentration.
        
        """

        self.NV = len(cobra_model.reactions)
        self.NEQ = len(reaction_indicies)
        self.reaction_indicies = np.asarray(reaction_indicies, dtype=np.int32)
        self.fva_tol = float(fva_tol)

        self._fva_dict = {'max' : {}, 'min' : {}}

        c_inner = np.array([r.objective_coefficient for r in
                            cobra_model.reactions])

        for index in self.reaction_indicies:

            # Initialize outer objectives
            c_maximize = np.zeros(self.NV)
            c_minimize = np.zeros(self.NV)
            c_maximize[index] = +1
            c_minimize[index] = -1

            # Initialize DFVA_core objects
            self._fva_dict['max'][index] = DFVA_Simulator_base(
                cobra_model, self.reaction_indicies, np.asarray(y0),
                bounds_func, c_inner, c_maximize, self.fva_tol)

            self._fva_dict['min'][index] = DFVA_Simulator_base(
                cobra_model, self.reaction_indicies, np.asarray(y0),
                bounds_func, c_inner, c_minimize, self.fva_tol)


    def integrate(self, t_start=0., t_end=1., res=200):
        """ Integrate the FVA problems from `t_start` to `t_end`, with `res`
        number of samples """

        self.res = res
        self.ts = np.linspace(t_start, t_end, res)
        
        for direction in ['max', 'min']:
            for index in self.reaction_indicies:
                self._fva_dict[direction][index].integrate(
                    t_start=float(t_start), t_end=float(t_end), res=int(res))
    
    def set_death_rate(self, death_rate):
        """ Set the death rate for each of the fva sub-problems """

        for direction in ['max', 'min']:
            for index in self.reaction_indicies:
                self._fva_dict[direction][index].set_death_rate(float(death_rate))
        
    def update_bounds_func_parameters(self, parameters):
        """ Update the parameters used by the attached ReactionBoundsFunction
        """

        for direction in ['max', 'min']:
            for index in self.reaction_indicies:
                self._fva_dict[direction][index].update_bounds_func_parameters(
                    np.asarray(parameters, dtype=np.float))

    def reset(self):
        """ Reset the integration to start at the y0 passed during class
        initialization """

        for direction in ['max', 'min']:
            for index in self.reaction_indicies:
                self._fva_dict[direction][index].reset()
    
    def _get_result(self, direction, result):

        property_shape = getattr(
            self._fva_dict[direction][self.reaction_indicies[0]], result).shape

        assert property_shape[1] == self.NEQ, "Result shape mismatch"

        out_arr = np.empty(property_shape)

        for i, index in enumerate(self.reaction_indicies):
            out_arr[:,i] = getattr(self._fva_dict[direction][index],
                                   result)[:,i]
        return out_arr

    @property
    def ys_max(self): return self._get_result('max', 'ys')

    @property
    def ys_min(self): return self._get_result('min', 'ys')

    @property
    def ys_dot_max(self): return self._get_result('max', 'ys_dot')

    @property
    def ys_dot_min(self): return self._get_result('min', 'ys_dot')

    # These are a bit more complicated: fluxes dont correspond directly to the
    # states that are tracked via FVA. This could find max/min of the vs fluxes
    # are are manipu
    # @property
    # def vs_max(self): return self._get_result('max', 'vs')

    # @property
    # def vs_min(self): return self._get_result('min', 'vs')




    # ys
    # ys_dot
    # vs


        






