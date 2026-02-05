from analysis_framework.Analysis import Analysis
from analysis_framework.Dataset import Dataset
from AltSetupHandler import AltSetupHandler
from itertools import combinations_with_replacement
import ROOT
import numpy as np
from pathlib import Path


class OptimalObservableHelper(Analysis):

    def __init__(self, dataset: Dataset, friend_datasets: list[Dataset] = []):
        super().__init__(dataset, friend_datasets)
        ROOT.gInterpreter.Declare("#include \"pol_helper.h\"")


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


    def define_optimal_observables_polarised(self, o_name: str, names: list[str], var_names: list[str], categories: list[str]|None = None, nom_epol: float = 0., nom_ppol: float = 0., add_pol_oo: bool = False):
        """Defines OOs, if names contains multiple names it averages over them."""
        nominal = [f"{name}_nominal" for name in names]
        nominal_pol_w = [f"(0.25 * (1 + {nom_epol} * gPolHelper.e_pol({hel_idx})) * (1 + {nom_ppol} * gPolHelper.p_pol({hel_idx})) * {n}[{hel_idx}])" for n in nominal for hel_idx in range(4)]
        combined_names = []
        for var_name in var_names:
            var = AltSetupHandler.get_var_from_name_1d(var_name)
            varied = [f"{name}_{var_name}" for name in names]
            # need to prefix every hel_sqme[hel_idx]
            # with 0.25 * (1 + nom_e_pol + gPolHelper.e_pol(hel_idx)) * (1 + nom_p_pol * gPolHelper.p_pol(hel_idx))
            # BUT I need to iterate hel_idx here myself and not use the global one because the OO is not allowed to know it!
            varied_pol_w = [f"(0.25 * (1 + {nom_epol} * gPolHelper.e_pol({hel_idx})) * (1 + {nom_ppol} * gPolHelper.p_pol({hel_idx})) * {v}[{hel_idx}])" for v in varied for hel_idx in range(4)]
            comb_name = f"{o_name}_{var_name}"
            self._define((comb_name, f"{1/var} * (({'+'.join(varied_pol_w)}) - ({'+'.join(nominal_pol_w)})) / ({'+'.join(nominal_pol_w)})"), categories)
            combined_names.append(comb_name)

        if add_pol_oo:
            # need to split nominal sqmes into helicities and average over the names
            sqme_hel_parts = []
            n = len(names)
            for i in range(4):
                sqme_hel_parts.append("+".join(nominal_pol_w[i*n:(i+1)*n]))
            self._define((f"{o_name}_epol", f"((({sqme_hel_parts[2]} + {sqme_hel_parts[3]}) / (1 + {nom_epol})) - (({sqme_hel_parts[0]} + {sqme_hel_parts[1]}) / (1 - {nom_epol}))) / ({'+'.join(sqme_hel_parts)})"), categories=categories)
            self._define((f"{o_name}_ppol", f"((({sqme_hel_parts[1]} + {sqme_hel_parts[3]}) / (1 + {nom_ppol})) - (({sqme_hel_parts[0]} + {sqme_hel_parts[2]}) / (1 - {nom_ppol}))) / ({'+'.join(sqme_hel_parts)})"), categories=categories)
            combined_names.append(f"{o_name}_epol")
            combined_names.append(f"{o_name}_ppol")
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


    def book_weight_sums(self, names: list[str], categories: list[str]|None = None, hel: bool = False) -> list[str]:
        if not hel:
            weight_names = [f"weight_{name}" for name in names]
        else:
            weight_names = [f"weight_hel_{name}" for name in names]
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


    @staticmethod
    def get_pol_var_list(weight_names: list[str], vary_par: str) -> list[float]:
        pol_var_list = []
        for w_name in weight_names:
            if vary_par in w_name:
                var_name = w_name.removeprefix("weight_").removeprefix("hel_")
                var = AltSetupHandler.get_var_from_name_1d(var_name)
                if var == 3e-4:
                    # FIXME: needs interface change :(
                    var = 2.5e-4
                elif var == -3e-4:
                    # FIXME: needs interface change :(
                    var = -2.5e-4
                pol_var_list.append(var)
        return pol_var_list


    def calc_oo_means(self, oo_names: list[str], weight_names: list[str], e_pol: float = 0.0, p_pol: float = 0.0, categories: list[str]|None = None, vary_pol: bool = False) -> dict[str, float]:
        oo_means = {}
        if vary_pol:
            # need to do very painful backward parsing because original interface did not have this in mind
            # hardcode this for now
            pol_vary_par = "g1z"
            pol_var_list = self.get_pol_var_list(weight_names, pol_vary_par)
        for w_name in weight_names:
            w_sum = self.get_sum(w_name, e_pol=e_pol, p_pol=p_pol, categories=categories)
            for oo in oo_names:
                oo_w_name = f"{oo}_{w_name}"
                oo_sum = self.get_sum(oo_w_name, e_pol=e_pol, p_pol=p_pol, categories=categories)
                print(oo,oo_w_name, oo_sum, w_sum)
                oo_means[oo_w_name] = oo_sum / w_sum #if w_sum != 0. else 0.

        if vary_pol:
            # print(pol_var_list)
            nom_w = [n for n in weight_names if "nominal" in n][0]
            for var in pol_var_list:
                w_sum_epol = self.get_sum(nom_w, e_pol=e_pol+var, p_pol=p_pol, categories=categories)
                w_sum_ppol = self.get_sum(nom_w, e_pol=e_pol, p_pol=p_pol+var, categories=categories)
                for oo in oo_names:
                    oo_w_name_epol = f"{oo}_weight_hel_{AltSetupHandler.make_name("epol", var)}"
                    oo_sum_epol = self.get_sum(f"{oo}_{nom_w}", e_pol=e_pol+var, p_pol=p_pol, categories=categories)
                    oo_means[oo_w_name_epol] = oo_sum_epol / w_sum_epol #if w_sum_epol != 0. else 0.

                    oo_w_name_ppol = f"{oo}_weight_hel_{AltSetupHandler.make_name("ppol", var)}"
                    oo_sum_ppol = self.get_sum(f"{oo}_{nom_w}", e_pol=e_pol, p_pol=p_pol+var, categories=categories)
                    oo_means[oo_w_name_ppol] = oo_sum_ppol / w_sum_ppol #if w_sum_ppol != 0. else 0.
        return oo_means


    @staticmethod
    def make_slope_graphs(oo_names: list[str], pars: list[str], x_points:
                          list[float], oo_means: dict[str, float], hel: bool = False) -> dict:
        graphs = {}
        w_name = "weight" if not hel else "weight_hel"
        for oo in oo_names:
            for par in pars:
                nominal = oo_means[f"{oo}_{w_name}_nominal"]
                x = x_points.copy()
                vars = [AltSetupHandler.make_name(par, p) for p in x_points]
                # y = [oo_means[f"{oo}_weight_{var}"] / nominal for var in vars]
                # y = [oo_means[f"{oo}_weight_{var}"] / nominal -1. for var in vars]
                y = [(oo_means[f"{oo}_{w_name}_{var}"] / nominal -1.) * 100  for var in vars]
                x.append(0.)
                # y.append(1.)
                y.append(0.)
                vals = sorted(zip(x, y), key=lambda a: a[0])
                # print(vals)
                g = ROOT.TGraph()
                for v in vals:
                    g.AddPoint(v[0], v[1])
                # g.SetTitle(f"{oo};{par}; #overline{{O}} / #overline{{O}}_{{0}}")
                graphs[f"{oo}_{par}"] = g
        return graphs


    @staticmethod
    def make_ratio_graph(graph, func):
        res_graph = ROOT.TGraph()
        for i in range(graph.GetN()):
            x = graph.GetPointX(i)
            y = graph.GetPointY(i)
            fx = func(x)
            # r_y = y-fx
            r_y = (y-fx) * 1e4 # only use 4 here as the values that get subtracted are already scaled by 100
            # r_y = (y / fx -1.) * 1e6
            # r_y = ((y+1.) / (fx+1.) -1.) * 1e6
            # r_y = (y / fx -1.) * 1e6 if fx != 0. else 0.
            res_graph.AddPoint(x, r_y)
        return res_graph


    @staticmethod
    def get_slopes(oo_names: list[str], pars: list[str], oo_means: dict[str, float], g: float = 1e-8, hel: bool = False) -> dict:
        slopes = {}
        w_name = "weight" if not hel else "weight_hel"
        for oo in oo_names:
            nominal = oo_means[f"{oo}_{w_name}_nominal"]
            for par in pars:
                var = AltSetupHandler.make_name(par, g)
                variation = oo_means[f"{oo}_{w_name}_{var}"]
                slope = (variation - nominal) / (g * nominal)
                slopes[f"{oo}_{par}"] = slope
        return slopes


    def print_fit_input(self, oo_names: list[str], ext_weight_names: list[str], pars: list[str], e_pol: float = 0.0, p_pol: float = 0.0, categories: list[str]|None = None, dir: str|None = None, name: str = "default", hel: bool =  False):
        # normalize everything to 1 ab_inv
        lumi = 1000
        w_name = "weight" if not hel else "weight_hel"
        n_events = self.get_sum(f"{w_name}_nominal", int_lumi=lumi, e_pol=e_pol, p_pol=p_pol, categories=categories)

        weight_names = ext_weight_names.copy()

        means = self.calc_oo_means(oo_names, weight_names, e_pol, p_pol, categories, vary_pol=hel)

        # TODO: also get other means here
        # need to find a way to somehow also add the pol variations to the weight_names
        pol_var_list = self.get_pol_var_list(weight_names, "g1z") if hel else []
        for var in pol_var_list:
            weight_names.append(f"weight_hel_{AltSetupHandler.make_name('epol', var)}")
            weight_names.append(f"weight_hel_{AltSetupHandler.make_name('ppol', var)}")

        means_dict = {w_name: [means[f"{oo}_{w_name}"] for oo in oo_names] for w_name in weight_names}


        # means_vec = np.asarray([means[f"{oo}_weight_nominal"] for oo in oo_names])
        means_vec = np.asarray(means_dict[f"{w_name}_nominal"])


        n_obs = len(oo_names)

        # TODO: replace with something more robust
        # pars = [oo.split("_")[-3] for oo in oo_names]
        n_pars = len(pars)
        slopes = self.get_slopes(oo_names, pars, means, hel=hel)
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
        text += "# pars, evt/ab_inv, means, slopes, cov\n"
        text += f"{pars}\n"
        text += f"{n_events}\n"
        text += f"{means_vec.tolist()}\n"
        text += f"{slope_mat.tolist()}\n"
        text += f"{C.tolist()}\n"
        text += f"{means_dict}\n"
        print(text)
        if dir:
            Path(dir).mkdir(parents=True, exist_ok=True)
            with open(f"{dir}/fit-inputs-{e_pol}-{p_pol}-{name}.txt", "w") as outfile:
                outfile.write(text)
