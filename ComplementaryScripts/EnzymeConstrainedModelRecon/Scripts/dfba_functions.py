#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/17/20

"""dfba_functions.py
:description : script
:param : 
:returns: 
:rtype:
TODO: unit?????
"""
import math

import cobra
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from tqdm import tqdm


class dfba_model(dict):

    def __init__(self, name):
        # dict.__init__(self)

        self['name'] = self.name = name
        self['cobra_model'] = None
        self['tspan'] = np.linspace(0, 10, 100)

        self['products_rxns'] = ['Biomass', ]
        self['initial_products'] = [0.1]

        self['substrate_rxns'] = ['EX_glc__D_e', ]
        self['initial_substrate'] = [10]
        self['substrate_bounds'] = ['lower_bound']

        self['opt_list'] = ['Biomass', ] + ['EX_glc__D_e', ]
        self['opt_list_direction'] = ['max'] * len(self['opt_list'])
        self['opt_list_multiply_by'] = np.array([1] * len(self['opt_list']))
        self['y0'] = None
        self['bound_method'] = 'actual_substrate_uptake'
        self['options'] = 'add_lp_feasibility'
        self['feasibility'] = 1000
        self['original_bounds'] = np.ones(len(self.substrate_rxns)) * 1000
        self['growth_tol'] = 1E-6

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(r"'Model' object has no attribute '%s'" % key)

    def __setattr__(self, key, value):
        if key in self.keys():
            self[key] = value
        else:
            raise AttributeError(r"'Model' object has no attribute '%s'" % key)

    def info(self):
        print('name: ', self.name, )
        print('bound_method: ', self.bound_method, )
        print('options: ', self.options, )
        print('opt_list: ', self.opt_list, )
        print('cobra_model: ', self.cobra_model, )

    def add_dynamic_bounds(self, y):
        """Use external concentrations to bound the uptake flux of glucose."""
        model = self.cobra_model
        biomass = y[0]

        if self.bound_method not in ['basic_Michaelis_Menten', 'actual_substrate_uptake']:
            print('dmodel.bound_method should be "basic_Michaelis_Menten" or "actual_substrate_uptake"')
            print('bound_method = "basic_Michaelis_Menten" setted')
            self.bound_method = 'basic_Michaelis_Menten'

        for index, rxn_id in enumerate(self.substrate_rxns):
            # substrate_co: substrate_concentration
            substrate_co = y[self.opt_list.index(rxn_id)]

            if self.bound_method == 'basic_Michaelis_Menten':  # TODO
                print('no code,please check')
                uptake_max_import = -10 * substrate_co / (5 + substrate_co)
                model.reactions.get_by_id(rxn_id).lower_bound = uptake_max_import

            elif self.bound_method == 'actual_substrate_uptake':  # actual_substrate_uptake = fba_sol flux
                # y_concentration = fluxes * biomass -->
                # fluxes(substrate_fluxes) = y_concentration/biomass
                substrate_fluxes = substrate_co / biomass  # * self.opt_list_multiply_by[index] # TODO：check lower_bound
                if self.substrate_bounds[index] == 'upper_bound':
                    uptake_max_import = min(substrate_fluxes, self.original_bounds[index])
                    model.reactions.get_by_id(rxn_id).upper_bound = uptake_max_import
                    # print(uptake_max_import)
                else:
                    uptake_max_import = max(-substrate_fluxes, self.original_bounds[index])
                    # uptake_max_import = -substrate_flux
                    model.reactions.get_by_id(rxn_id).lower_bound = uptake_max_import

    def rest(self):
        self.feasibility = 1000  # reste feasibility
        original_bounds = self.original_bounds
        for index, rxn_id in enumerate(self.substrate_rxns):
            if self.substrate_bounds[index] == 'upper_bound':
                original_bound = self.cobra_model.reactions.get_by_id(rxn_id).upper_bound
            else:
                original_bound = self.cobra_model.reactions.get_by_id(rxn_id).lower_bound
            original_bounds[index] = original_bound
        original_bounds[original_bounds > 1000] = 1000
        print('original_bounds:', original_bounds)
        self.original_bounds = original_bounds


def infeasible_event(t, y, d_model, epsilon=1E-6):
    """
    Determine solution feasibility.

    Avoiding infeasible solutions is handled by solve_ivp's built-in event detection.
    This function re-solves the LP to determine whether or not the solution is feasible
    (and if not, how far it is from feasibility). When the sign of this function changes
    from -epsilon to positive, we know the solution is no longer feasible.

    """

    if d_model.options == 'add_lp_feasibility':
        model = d_model.cobra_model
        with model:
            d_model.add_dynamic_bounds(y)

            cobra.util.add_lp_feasibility(model)
            feasibility = cobra.util.fix_objective_as_constraint(model)
        return feasibility - epsilon
    else:

        return d_model.feasibility


def dy_dt_system(t, y, d_model):
    """Calculate the time derivative of external species.
        y: concentration
    """
    biomass = y[0]

    model = d_model.cobra_model
    # print(model.objective.expression)
    # model.objective = d_model.opt_list[0]
    # Calculate the specific exchanges fluxes at the given external concentrations.
    with model:

        d_model.add_dynamic_bounds(y)
        # if d_model.options:
        if d_model.options == 'add_lp_feasibility':
            cobra.util.add_lp_feasibility(model)
            feasibility = cobra.util.fix_objective_as_constraint(model)
            lex_constraints = cobra.util.add_lexicographic_constraints(
                model, d_model.opt_list, d_model.opt_list_direction)
            fluxes = lex_constraints.values


        else:
            try:
                if d_model.options == 'only_FBA':
                    sol_fba = model.optimize()
                    fluxes = sol_fba.fluxes[d_model.opt_list].values
                else:
                    # feasibility = cobra.util.fix_objective_as_constraint(model)
                    lex_constraints = cobra.util.add_lexicographic_constraints(
                        model, d_model.opt_list, d_model.opt_list_direction)
                    fluxes = lex_constraints.values

                    # print(fluxes)
                d_model.feasibility = fluxes[0]
                # print(fluxes)
            except:
                fluxes = np.zeros(len(d_model.opt_list))
                # print('model optimize failed ')
                d_model.feasibility = 0
            if model.solver.status == 'infeasible' or d_model.feasibility <= d_model.growth_tol * 1e-3:
                d_model.feasibility = 0

    # Since the calculated fluxes are specific rates, we multiply them by the
    # biomass concentration to get the bulk exchange rates.

    y_concentration = fluxes * biomass * d_model.opt_list_multiply_by

    # This implementation is **not** efficient, so I display the current
    # simulation time using a progress bar.
    if dy_dt_system.pbar is not None:
        dy_dt_system.pbar.update(1)
        dy_dt_system.pbar.set_description('t = {:.3f}'.format(t))

    return y_concentration


def do_dfba_solve_ivp(d_model, reset=True, bool_tqdm=True, terminal=True):
    print('\n---------- Dynamic Flux Balance Analysis (dFBA) solve_ivp ... ---------- ')

    if reset:
        d_model.rest()

    dy_dt_system.pbar = None
    epsilon = d_model.growth_tol
    infeasible_event2 = lambda t, y: infeasible_event(t, y, d_model, epsilon=epsilon)
    infeasible_event2.epsilon = epsilon
    infeasible_event2.direction = 1
    infeasible_event2.terminal = terminal

    if bool_tqdm:
        with tqdm() as pbar:
            dy_dt_system.pbar = pbar

            sol = solve_ivp(
                fun=lambda t, y: dy_dt_system(t, y, d_model),
                events=infeasible_event2,
                t_span=(d_model.tspan.min(), d_model.tspan.max()),
                y0=d_model.y0,
                t_eval=d_model.tspan,
                # dense_output=True,
                rtol=1e-3,
                atol=1e-6,
                method='BDF')
    else:
        sol = solve_ivp(
            fun=lambda t, y: dy_dt_system(t, y, d_model),
            events=infeasible_event2,
            t_span=(d_model.tspan.min(), d_model.tspan.max()),
            y0=d_model.y0,
            t_eval=d_model.tspan,
            # dense_output=True,
            rtol=1e-3,
            atol=1e-6,
            method='BDF')

    result_sol = np.concatenate((sol.t.reshape((-1, 1)), sol.y.T), axis=1)
    return result_sol


def do_dfba_odeint(d_model, reset=True, bool_tqdm=True, ):
    print('\n---------- Dynamic Flux Balance Analysis (dFBA) odeint ... ---------- ')

    if reset:
        d_model.rest()

    dy_dt_system.pbar = None
    if bool_tqdm:
        with tqdm() as pbar:
            dy_dt_system.pbar = pbar

            sol = odeint(
                func=dy_dt_system,
                tfirst=True,
                y0=d_model.y0,
                t=d_model.tspan,
                args=(d_model,))
    else:
        sol = odeint(
            func=dy_dt_system,
            tfirst=True,
            y0=d_model.y0,
            t=d_model.tspan,
            args=(d_model,))

    result_sol = np.concatenate((d_model.tspan.reshape((-1, 1)), sol), axis=1)
    return result_sol


def do_dfba_iteration(d_model, reset=True, bool_tqdm=True):
    print('\n---------- Dynamic Flux Balance Analysis (dFBA) iteration ... ---------- ')

    if reset:
        d_model.rest()
    model = d_model.cobra_model

    # NOTE: ignore check part and other exchange reactions

    biomass = d_model.initial_products[0]
    t_step = (d_model.tspan.max() - d_model.tspan.min()) / (len(d_model.tspan) - 1)

    # Initialize bounds
    uptake_bound = d_model.initial_substrate / (biomass * t_step);

    # uptakeBound[uptakeBound>1000] = 1000
    original_bounds = d_model.original_bounds
    for index, rxn_id in enumerate(d_model.substrate_rxns):
        if d_model.substrate_bounds[index] == 'upper_bound':
            model.reactions.get_by_id(rxn_id).upper_bound = min(uptake_bound[index], original_bounds[index])
        else:
            model.reactions.get_by_id(rxn_id).lower_bound = max(-uptake_bound[index], original_bounds[index])

    # sol_products = [d_model.initial_products]
    # sol_substrate = []
    sol_opt_list = np.array(np.array(d_model.y0), ).reshape(1, -1)
    tspan = d_model.tspan

    with tqdm() as pbar:
        for t in tspan:

            try:
                sol = model.optimize()
            except:
                print('\nNo feasible solution - nutrients exhausted. Biomass:\t %f\n', biomass)
                break

            mu = sol.objective_value

            if mu <= d_model.growth_tol * (1e-3):
                print('\nNo feasible solution - nutrients exhausted. Biomass:\t %f\n', biomass)
                break
            # uptakeFlux = sol.fluxes[d_model.substrate_rxns].values

            y = sol.fluxes[d_model.opt_list].values
            y = y * d_model.opt_list_multiply_by
            biomass = biomass * math.exp(mu * t_step)
            y[0] = biomass
            y[1:] = sol_opt_list[-1, 1:] - y[1:] / mu * biomass * (1 - math.exp(mu * t_step))

            if y.min() < 0:
                break

            # sol_opt_list = np.append(sol_opt_list,y, axis=0)
            sol_opt_list = np.vstack([sol_opt_list, y])

            # Update bounds for uptake reactions
            uptake_bound = y / (biomass * t_step)
            # uptake_bound = sol.fluxes[d_model.substrate_rxns].values / (biomass * t_step)

            # This is to avoid any numerical issues
            uptake_bound[uptake_bound > 1000] = 1000
            # Figure out if the computed bounds were above the original bounds
            # aboveOriginal = (uptakeBound > originalBound) & (originalBound > 0);
            # Revert to original bounds if the rate was too high
            # uptakeBound(aboveOriginal) = originalBound(aboveOriginal);
            uptake_bound[abs(uptake_bound) <= d_model.growth_tol] = 0
            for index, rxn_id in enumerate(d_model.substrate_rxns):
                # index = d_model.substrate_rxns.index(rxn_id)
                if d_model.substrate_bounds[index] == 'upper_bound':
                    model.reactions.get_by_id(rxn_id).upper_bound = min(uptake_bound[index], original_bounds[index])
                    # print(min(uptake_bound[index], original_bounds[index]))
                else:
                    model.reactions.get_by_id(rxn_id).lower_bound = max(-uptake_bound[index], original_bounds[index])

            pbar.update(1)
            pbar.set_description('t = {:.3f}'.format(t))

    actual_tspan = tspan[tspan <= t]
    result_sol = np.concatenate((actual_tspan.reshape((-1, 1)), sol_opt_list[0:actual_tspan.shape[0], :]), axis=1)

    return result_sol


def draw(result_sol, d_model):
    fig, axs = plt.subplots(1, 1)
    markersize = 4
    ax = axs
    line_1 = ax.plot(result_sol[:, 0], result_sol[:, 1], '^-', markersize=markersize, color='r', alpha=0.5)

    ax2 = plt.twinx(ax)
    for i in range(2, result_sol.shape[1]):
        line_2 = ax2.plot(result_sol[:, 0], result_sol[:, i], '.-', markersize=markersize, alpha=0.5)

    ax.set_xlim(min(d_model.tspan), max(d_model.tspan))
    ax.set_ylabel('Biomass ', color='r')
    ax2.set_ylabel('Metabolites ', color='black')
    ax.legend(['biomass'], loc=2)
    ax2.legend(d_model.opt_list[1:], loc=1)
    plt.show()
    plt.close()


def get_batch_experiment_data(draw=True):
    df = pd.read_csv('../Fermentation/batch.csv')  # mmol/L  = mM
    df['PHA (mM)'] = df['PHA'] * 1000 / 150  # TODO: check
    df['citrate (mM)'] = df['citrate (g/L)'] * 1000 / 189
    df['glucose (mM)'] = df['glucose(g/L)'] * 1000 / 180
    df['pyruvate (mM)'] = df['pyruvate(g/L)'] * 1000 / 87
    df['acecate (mM)'] = df['acecate(g/L)'] * 1000 / 59

    mets = ['PHA (mM)', 'citrate (mM)', 'glucose (mM)', 'pyruvate (mM)', 'acecate (mM)']
    if draw:
        fig, axs = plt.subplots(2, 1)
        markersize = 4
        ax = axs[0]
        ax2 = axs[1]  # plt.twinx(ax)
        ax.plot(df['time(h)'], df['X(gdw/L)'], '--o', markersize=markersize)
        for i in mets:
            ax2.plot(df['time(h)'], df[i], '--o', markersize=markersize, )

        ax.set_ylabel('Biomass')
        ax2.set_ylabel('Metabolites (mM)')
        ax.legend(['biomass'], loc=0)
        ax2.legend(mets, loc=0)
        plt.show()
        plt.close()

    return df


if __name__ == '__main__':
    # %% <local test>
    case = 3
    if case == 1:
        from cobra.test import create_test_model

        model = create_test_model('textbook')
        biomass_id = 'Biomass_Ecoli_core'
    elif case == 2:
        model = cobra.io.read_sbml_model('../../../ComplementaryData/iML1515.xml')
        biomass_id = 'BIOMASS_Ec_iML1515_core_75p37M'
    if case == 3:
        model = cobra.io.load_json_model('test_dfba_model.json')
        biomass_id = 'Biomass_v1'

    d_model = dfba_model('test_model')
    d_model.cobra_model = model.copy()
    d_model.tspan = np.linspace(0, 10, 101)
    d_model.products_rxns = [biomass_id, 'EX_co2_e', 'EX_ac_e', 'EX_pyr_e']
    # d_model.initial_products = [0.1, 0]
    d_model.substrate_rxns = ['EX_glc__D_e']
    # d_model.initial_substrate = [10]
    d_model.opt_list = [biomass_id, 'EX_co2_e', 'EX_ac_e', 'EX_pyr_e', 'EX_glc__D_e']
    d_model.y0 = [0.1, 0, 0, 0, 10]
    d_model.opt_list_direction = ['max', ] * len(d_model.opt_list)
    d_model.opt_list_multiply_by = np.array([1] * len(d_model.opt_list))
    d_model.options = 'only_FBA'  # if add_lp_feasibility(default) will be slow
    d_model.bound_method = 'actual_substrate_uptake'  # 'basic_Michaelis_Menten'
    d_model.growth_tol = 1E-6

    if case == 3:
        d_model.tspan = np.linspace(0, 50, 101)

        glc_id = 'Exchange_Glucopyranose'
        d_model.products_rxns = [biomass_id, 'PHA_secretion',
                                 'Exchange_ACET']  # 'PHA_secretion', 'Exchange_ACET', 'Exchange_PYRUVATE',
        # 'PHA_secretion' 'PHA_secretion'
        d_model.substrate_rxns = [glc_id]
        d_model.substrate_bounds = ['upper_bound']
        d_model.opt_list = d_model.products_rxns + d_model.substrate_rxns  # [biomass_id, 'EX_co2_e', 'EX_ac_e', 'EX_pyr_e', 'EX_glc__D_e']
        d_model.opt_list_direction = ['max', ] * (len(d_model.opt_list) - 1) + ['min']
        d_model.opt_list_multiply_by = np.array([1] * (len(d_model.opt_list) - 1) + [-1])
        d_model.y0 = [0.07446502] + [0] * (len(d_model.opt_list) - 2) + [28]
        d_model.options = 'not_add_lp_feasibility'  # not_add_lp_feasibility
        d_model.growth_tol = 1E-6

    # %%
    d_model.cobra_model = model.copy()
    sol_iteration = do_dfba_iteration(d_model)
    draw(sol_iteration, d_model)

    # %% <do_dfba_solve_ivp>
    d_model.cobra_model = model.copy()
    sol_solve_ivp = do_dfba_solve_ivp(d_model, bool_tqdm=False, terminal=True)
    draw(sol_solve_ivp, d_model)
    # %% <do_dfba_odeint>
    d_model.cobra_model = model.copy()
    sol_odeint = do_dfba_odeint(d_model, bool_tqdm=False, )
    draw(sol_odeint, d_model)
