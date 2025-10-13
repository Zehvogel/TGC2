from analysis_framework.Analysis import Analysis
from AltSetupHandler import AltSetupHandler
from itertools import combinations_with_replacement
import ROOT
import numpy as np
from pathlib import Path


class OptimalObservableHelper(Analysis):

    def __init__(self, dataset):
        super().__init__(dataset)


    def define_optimal_observables(self, o_name: str, names: list[str], var_names: list[str], categories: list[str]|None = None):
        """Defines OOs, if names contains multiple names it averages over them."""
        nominal = [f"{name}_nominal" for name in names]
        combined_names = []
        for var_name in var_names:
            var = AltSetupHandler.get_var_from_name_1d(var_name)
            varied = [f"{name}_{var_name}" for name in names]
            comb_name = f"{o_name}_{var_name}"
            self._define((comb_name, f"{1/var} * (({'+'.join(varied)}) - ({'+'.join(nominal)})) / ({'+'.join(nominal)})"), categories)
            combined_names.append(comb_name)
        return combined_names


    # exists both here and in the reweight helper :/
    def define_optimal_observables_truth(self, names: list[str], truth_categories: list[str]):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self._define((f"mc_O_{name}", f"{1/var} * (mc_sqme_{name} - mc_sqme_nominal) / mc_sqme_nominal"), truth_categories)


    def define_optimal_observables_truth_averaged(self, names: list[str], truth_categories: list[str]):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self._define((f"av_mc_O_{name}", f"{1/var} * ((mc_sqme_{name} + wj_mc_sqme_{name}) - (mc_sqme_nominal + wj_mc_sqme_nominal)) / (mc_sqme_nominal + wj_mc_sqme_nominal)"), truth_categories)


    def define_optimal_observables_reco(self, names: list[str]):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self.Define(f"O_{name}", f"{1/var} * ((reco_sqme_12_{name} + reco_sqme_21_{name}) - (reco_sqme_12_nominal + reco_sqme_21_nominal)) / (reco_sqme_12_nominal + reco_sqme_21_nominal)")


    def book_oo_matrix(self, observables: list[str], categories: list[str]|None = None):
        for o1, o2 in combinations_with_replacement(observables, 2):
            name = f"{o1}_{o2}"
            self._define((name, f"{o1} * {o2}"), categories=categories)
            self.book_sum(name, name, categories=categories)


    def get_oo_matrix(self, observables: list[str], int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0, categories: list[str]|None = None):
        res = []
        for o1, o2 in combinations_with_replacement(observables, 2):
            name = f"{o1}_{o2}"
            res.append(self.get_sum(name, int_lumi=int_lumi, e_pol=e_pol, p_pol=p_pol, categories=categories))
        return res


    def book_weight_sums(self, names: list[str], categories: list[str]|None = None) -> list[str]:
        weight_names = [f"weight_{name}" for name in names]
        for w_name in weight_names:
            self.book_sum(w_name, w_name, categories=categories)
        return weight_names


    def define_weighted_oo(self, oo_names: list[str], weight_names: list[str], categories: list[str]|None = None):
        for oo in oo_names:
            for w_name in weight_names:
                oo_w_name = f"{oo}_{w_name}"
                self._define((oo_w_name, f"{oo} * {w_name}"), categories=categories)


    def book_oo_sums(self, oo_names: list[str], weight_names: list[str], categories: list[str]|None = None):
        for oo in oo_names:
            for w_name in weight_names:
                oo_w_name = f"{oo}_{w_name}"
                self.book_sum(oo_w_name, oo_w_name, categories=categories)


    def book_oo_histograms(self, oo_names: list[str], weight_names: list[str], categories: list[str]|None = None):
        for oo in oo_names:
            for w_name in weight_names:
                oo_w_name = f"{oo}_{w_name}"
                self.book_histogram_1D(oo_w_name, oo_w_name, ("", "", 250, -5, 5), categories=categories)


    def calc_oo_means(self, oo_names: list[str], weight_names: list[str], e_pol: float = 0.0, p_pol: float = 0.0, categories: list[str]|None = None) -> dict[str, float]:
        oo_means = {}
        for w_name in weight_names:
            w_sum = self.get_sum(w_name, e_pol=e_pol, p_pol=p_pol, categories=categories)
            for oo in oo_names:
                oo_w_name = f"{oo}_{w_name}"
                oo_sum = self.get_sum(oo_w_name, e_pol=e_pol, p_pol=p_pol, categories=categories)
                oo_means[oo_w_name] = oo_sum / w_sum
        return oo_means


    def make_slope_graphs(self, oo_names: list[str], pars: list[str], alt_config_names:
                          list[str], oo_means: dict[str, float]) -> dict:
        graphs = {}
        for oo in oo_names:
            for par in pars:
                vars = [var for var in alt_config_names if par in var]
                nominal = oo_means[f"{oo}_weight_nominal"]
                x = [AltSetupHandler.get_var_from_name_1d(var) for var in vars]
                y = [oo_means[f"{oo}_weight_{var}"] / nominal for var in vars]
                x.append(0.)
                y.append(1.)
                vals = sorted(zip(x, y), key=lambda a: a[0])
                # print(vals)
                g = ROOT.TGraph()
                for v in vals:
                    g.AddPoint(v[0], v[1])
                g.SetTitle(f"{oo};{par}; #overline{{O}} / #overline{{O}}_{{0}}")
                graphs[f"mc_{oo}_{par}"] = g
        return graphs

    @staticmethod
    def get_slopes(oo_names: list[str], pars: list[str], oo_means: dict[str, float]) -> dict:
        slopes = {}
        for oo in oo_names:
            nominal = oo_means[f"{oo}_weight_nominal"]
            for par in pars:
                var = f"{par}_pos_1em08"
                variation = oo_means[f"{oo}_weight_{var}"]
                g = 1e-8
                slope = (variation - nominal) / (g * nominal)
                slopes[f"{oo}_{par}"] = slope
        return slopes


    def print_fit_input(self, oo_names: list[str], weight_names: list[str], e_pol: float = 0.0, p_pol: float = 0.0, categories: list[str]|None = None, dir: str|None = None, name: str = "default"):
        # normalize everything to 1 ab_inv
        lumi = 1000
        n_events = self.get_sum("weight_nominal", int_lumi=lumi, e_pol=e_pol, p_pol=p_pol, categories=categories)

        means = self.calc_oo_means(oo_names, weight_names, e_pol, p_pol, categories)
        means_vec = np.asarray([means[f"{oo}_weight_nominal"] for oo in oo_names])

        n_obs = len(oo_names)

        # TODO: replace with something more robust
        pars = [oo.split("_")[-3] for oo in oo_names]
        n_pars = len(pars)
        slopes = self.get_slopes(oo_names, pars, means)
        slope_list = list(slopes.values())
        slope_mat = np.zeros((n_obs, n_pars))
        for i in range(n_obs):
            for j in range(n_pars):
                slope_mat[i, j] = slope_list[i*n_pars + j]

        mat = self.get_oo_matrix(oo_names, int_lumi=lumi, e_pol=e_pol, p_pol=p_pol, categories=categories)
        C_tilde = np.zeros((n_obs, n_obs))
        for k, (i, j) in enumerate(combinations_with_replacement(range(n_obs), 2)):
            C_tilde[i, j] = mat[k]
            C_tilde[j, i] = mat[k]

        C = C_tilde / n_events - np.outer(means_vec, means_vec)

        text = ""
        text += "#evt/ab_inv, means, slopes, cov\n"
        text += f"{n_events}\n"
        text += f"{means_vec.tolist()}\n"
        text += f"{slope_mat.tolist()}\n"
        text += f"{C.tolist()}\n"
        print(text)
        if dir:
            Path(dir).mkdir(parents=True, exist_ok=True)
            with open(f"{dir}/fit-inputs-{e_pol}-{p_pol}-{name}.txt", "w") as outfile:
                outfile.write(text)
