import ROOT
from ast import literal_eval
from itertools import combinations_with_replacement
import numpy as np
from AltSetupHandler import AltSetupHandler
from pathlib import Path


def load_fit_input(input_path):
    with open(input_path) as infile:
        lines = infile.readlines()
    pars = literal_eval(lines[1])
    n_per_ab = literal_eval(lines[2])
    nominal_means = literal_eval(lines[3])
    slopes = literal_eval(lines[4])
    Cov_O = literal_eval(lines[5])
    means_list = literal_eval(lines[6])
    # print(lines)
    return pars, n_per_ab, nominal_means, slopes, Cov_O, means_list


def build_model(workspace, parameters, obs, nominal_means, slopes, Cov_mean_O, rename_var_suffix: str = ""):
    mus_exp = []
    all_vars = []
    rel_changes = []
    for i, m in enumerate(nominal_means):
        vars = []
        for j in range(len(parameters)):
            var = ROOT.RooProduct(f"var_{i}_{j}{rename_var_suffix}", f"var_{i}_{j}{rename_var_suffix}", parameters[j], ROOT.RooFit.RooConst(slopes[i][j]))
            vars.append(var)
        all_vars += vars
        nominal_exp = ROOT.RooFit.RooConst(m)
        # mu_exp = ROOT.RooAddition(f"mu_exp_{i}", f"mu_exp_{i}", [nominal_exp] + vars)
        rel_change = ROOT.RooAddition(f"rel_change_{i}{rename_var_suffix}", f"rel_change_{i}{rename_var_suffix}", [ROOT.RooFit.RooConst(1.0)] + vars)
        rel_changes.append(rel_change)
        mu_exp = ROOT.RooProduct(f"mu_exp_{i}{rename_var_suffix}", f"mu_exp_{i}{rename_var_suffix}", rel_change, nominal_exp)
        mus_exp.append(mu_exp)

    # print(all_vars)

    n_obs = len(nominal_means)
    Cov_mean_O_root = ROOT.TMatrixDSym(n_obs)
    for i, j in combinations_with_replacement(range(n_obs), r=2):
        Cov_mean_O_root[i][j] = Cov_mean_O[i, j]
        Cov_mean_O_root[j][i] = Cov_mean_O[i, j]

    fit_model = ROOT.RooMultiVarGaussian(f"fit_model{rename_var_suffix}", f"fit_model{rename_var_suffix}", obs, mus_exp, Cov_mean_O_root)
    workspace.Import(fit_model, Silence=True)


def create_models(workspace, pars, nominal_means, slopes, Cov_mean_O):
    parameters = [ROOT.RooRealVar(par, par, 0., -0.5, 0.5) for par in pars]

    obs = []
    for i in range(len(nominal_means)):
        sigma = np.sqrt(Cov_mean_O[i,i])
        # only really meaningful for the Asimov case where obs = means, for the datasets where the parameters vary they could easily be outside this range if stats are large
        # is it still valid then to use the constant covariance matrix without scaling it itself? The fit still seems to work :)
        o = ROOT.RooRealVar(f"O_{pars[i]}", f"O_{pars[i]}", nominal_means[i], nominal_means[i]-500*sigma, nominal_means[i]+500*sigma)
        obs.append(o)

    build_model(workspace, parameters, obs, nominal_means, slopes, Cov_mean_O, rename_var_suffix="")

    # mean_vec = ROOT.TVectorD(n_obs)
    # for i in range(n_obs):
    #     mean_vec[i] = nominal_means[i]

    # gen_model = ROOT.RooMultiVarGaussian("gen_model", "gen_model", obs, mean_vec, Cov_mean_O_root)
    # workspace.Import(gen_model)

    workspace.saveSnapshot("nominal_parameters", parameters)
    workspace.saveSnapshot("nominal_observables", obs)
    workspace.defineSet("observables", obs)
    workspace.defineSet("parameters", parameters)


def build_sim_ws(runs, input_path, oo_name, run_ab, split_helicity_reversal: bool = True):
    w = ROOT.RooWorkspace("w")
    pars_list = []
    n_per_ab = []
    nominal_means = []
    slopes = []
    Cov_mean_O = []
    means_list = []
    for lumi_share, e_pol, p_pol in runs:
        fit_input_file = f"{input_path}/fit-inputs-{e_pol}-{p_pol}-{oo_name}_oo.txt"
        fit_inputs = load_fit_input(fit_input_file)
        pars_list.append(fit_inputs[0])
        n_per_ab.append(fit_inputs[1])
        nominal_means.append(fit_inputs[2])
        slopes.append(fit_inputs[3])
        Cov_O = fit_inputs[4]
        Cov_mean_O.append(np.asarray(Cov_O) / (n_per_ab[-1] * run_ab * lumi_share))
        means_list.append(fit_inputs[5])

    # build common obs and pars, assume they are the same for all runs
    n = len(nominal_means[0])
    # max_sigmas = [0.] * n
    # for i in range(n):
    #     for j in range(len(runs)):
    #         sigma = np.sqrt(Cov_mean_O[j][i,i])
    #         if sigma > max_sigmas[i]:
    #             max_sigmas[i] = sigma
    obs = []
    for i in range(n):
        # o = ROOT.RooRealVar(f"O_{pars_list[0][i]}", f"O_{pars_list[0][i]}", nominal_means[0][i], nominal_means[0][i]-5*max_sigmas[i], nominal_means[0][i]+5*max_sigmas[i])
        # super annoying to do correctly so just choose something big
        o = ROOT.RooRealVar(f"O_{pars_list[0][i]}", f"O_{pars_list[0][i]}", nominal_means[0][i], -5., 5.)
        obs.append(o)
    pars = pars_list[0].copy()
    if split_helicity_reversal:
        pars = pars[:-2]  # remove pol pars
    fit_parameters = [ROOT.RooRealVar(par, par, 0., -0.5, 0.5) for par in pars]
    e_pol_fit_parameters = []
    p_pol_fit_parameters = []
    if split_helicity_reversal:
        e_pol_fit_parameters = [
            ROOT.RooRealVar("e_pol_L", "e_pol_L", 0., -0.5, 0.5),
            ROOT.RooRealVar("e_pol_R", "e_pol_R", 0., -0.5, 0.5)
            ]
        p_pol_fit_parameters = [
            ROOT.RooRealVar("p_pol_L", "p_pol_L", 0., -0.5, 0.5),
            ROOT.RooRealVar("p_pol_R", "p_pol_R", 0., -0.5, 0.5)
            ]
    all_fit_parameters = fit_parameters + e_pol_fit_parameters + p_pol_fit_parameters


    models = []
    for i, r_conf in enumerate(runs):
        lumi_share, e_pol, p_pol = r_conf
        nominal_means_run = nominal_means[i]
        slopes_run = slopes[i]
        Cov_mean_O_run = Cov_mean_O[i]
        # adjust slopes for pol parameters if split
        run_fit_parameters = fit_parameters.copy()
        if split_helicity_reversal:
            if e_pol > 0:
                e_pol_par = e_pol_fit_parameters[1]  # R
            else:
                e_pol_par = e_pol_fit_parameters[0]  # L
            if p_pol > 0:
                p_pol_par = p_pol_fit_parameters[1]  # R
            else:
                p_pol_par = p_pol_fit_parameters[0]  # L
            run_fit_parameters += [e_pol_par, p_pol_par]
        build_model(w, run_fit_parameters, obs, nominal_means_run, slopes_run, Cov_mean_O_run, f"_run{i}")
        models.append(w.pdf(f"fit_model_run{i}"))

    run_categories = ROOT.RooCategory("runs", "Run configurations", {f"run_{i}": i for i in range(len(models))})
    model_map = {f"run_{i}": models[i] for i in range(len(models))}
    sim_model = ROOT.RooSimultaneous("sim_model", "sim_model", model_map, run_categories)
    w.Import(sim_model)
    w.defineSet("observables", obs)
    w.saveSnapshot("nominal_parameters", all_fit_parameters)
    return w



def make_datasets(workspace, prefix: str, pars: list[str], means_dict, x_points: list[float] = []):
    obs = workspace.set("observables")
    res = []
    names = [f"{prefix}_nominal"]
    names += [f"{prefix}_{AltSetupHandler.make_name(par, x)}" for par in pars for x in x_points]
    print(names)
    for name in names:
        means = means_dict[name]
        for i, p in enumerate(pars):
            obs[f"O_{p}"] = means[i]
            # print(f"set O_{p} to {means[i]}")
        ds = ROOT.RooDataSet(name, name, obs)
        ds.add(obs)
        # print(ds)
        # workspace.Import(ds)
        # FIXME: workspace bug can not store all datasets
        # https://github.com/root-project/root/issues/20904
        res.append(ds)
    return res


# def fit(workspace, model_name, ds_name):
#     workspace.loadSnapshot("nominal_parameters")
#     model = workspace.pdf(model_name)
#     ds = workspace.data(ds_name)
#     model.fitTo(ds)
#     params = workspace.set("parameters")
#     for p in params:
#         print(p.GetName(), p.getVal(), p.getError())


def fit(workspace, model_name, ds, silent: bool = False, minos: bool = False):
    """Fit the model to the given dataset and return the fitted parameters and their errors."""
    workspace.loadSnapshot("nominal_parameters")
    model = workspace.pdf(model_name)
    # could also store fit result here
    # model.fitTo(ds)
    # res = model.fitTo(ds, Minos=True, Save=True, PrintLevel=-1 if silent else 1)
    # minos is not really needed here and takes a lot of time
    res = model.fitTo(ds, Minos=minos, Save=True, PrintLevel=-1 if silent else 1)
    return res


def make_par_histo(fit_res, parameter_names: list[str]):
    n = len(parameter_names)
    h = ROOT.TH1D("", "", n, 0, n)
    fit_parameters = fit_res.floatParsFinal()
    for name in parameter_names:
        try:
            par = fit_parameters.find(name)
            h.Fill(name, par.getError())
        except:
            h.Fill(name, 0)
    return h


def make_par_histos(fit_results, histos, stack, l, parameter_names: list[str], legend_name_dict: dict[str, str] = {}):
    for i, (name, fit_res) in enumerate(fit_results.items()):
        h = make_par_histo(fit_res, parameter_names)
        h.SetFillColor(ROOT.kP10Blue + i)
        leg_name = legend_name_dict[name] if name in legend_name_dict else name
        l.AddEntry(h, leg_name, "f")
        histos[name] = h
        stack.Add(h)


def make_plots_runs_per_oo_name(oo_name: str, fit_results: dict,
                                parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"],
                                legend_pars = None, legend_name_dict: dict[str, str] = {}):
    histos = {}
    stack = ROOT.THStack()
    if legend_pars:
        l = ROOT.TLegend(*legend_pars)
    else:
        l = ROOT.TLegend()
    make_par_histos(fit_results[oo_name], histos, stack, l, parameter_names, legend_name_dict)
    return histos, stack, l


def make_plots_oo_names_per_run_name(run_name: str, fit_results: dict, oo_names: list[str], parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"], legend_pars = None, legend_name_dict: dict[str, str] = {}):
    histos = {}
    stack = ROOT.THStack()
    if legend_pars:
        l = ROOT.TLegend(*legend_pars)
    else:
        l = ROOT.TLegend()
    for i, oo_name in enumerate(oo_names):
        fit_res = fit_results[oo_name][run_name]
        h = make_par_histo(fit_res, parameter_names)
        h.SetFillColor(ROOT.kP10Blue + i)
        leg_name = legend_name_dict[oo_name] if oo_name in legend_name_dict else oo_name
        l.AddEntry(h, leg_name, "f")
        histos[oo_name] = h
        stack.Add(h)
    return histos, stack, l


class Plotter:
    def __init__(self):
        self.histos = {}
        self.stacks = {}
        self.legends = {}
        self.canvases = {}


    def apply_canvas_properties(self, c):
        c.SetLeftMargin(0.15)
        c.SetRightMargin(0.05)
        c.SetBottomMargin(0.1)


    def apply_stack_properties(self, stack):
        # stack.GetXaxis().SetTitleSize(0.5)
        stack.GetXaxis().SetLabelSize(0.125)
        stack.GetYaxis().SetTitle("absolute uncertainty")


    def draw_plots_runs_per_oo_name(self, plot_name: str, oo_name: str, fit_results: dict,
                                parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"],
                                legend_pars = None, legend_name_dict: dict[str, str] = {}):
        histos, stack, l = make_plots_runs_per_oo_name(oo_name, fit_results, parameter_names, legend_pars, legend_name_dict)
        self.histos[plot_name] = histos
        self.stacks[plot_name] = stack
        self.legends[plot_name] = l
        c = ROOT.TCanvas()
        stack.Draw("nostackb hist")
        self.apply_stack_properties(stack)
        l.Draw()
        self.apply_canvas_properties(c)
        c.Draw()
        self.canvases[plot_name] = c


    def draw_plots_oo_names_per_run_name(self, plot_name: str, run_name: str, fit_results: dict, oo_names: list[str],
                                parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"],
                                legend_pars = None, legend_name_dict: dict[str, str] = {}, plot_dir: str|None = None):
        histos, stack, l = make_plots_oo_names_per_run_name(run_name, fit_results, oo_names, parameter_names, legend_pars, legend_name_dict)
        self.histos[plot_name] = histos
        self.stacks[plot_name] = stack
        self.legends[plot_name] = l
        if plot_dir:
            Path(plot_dir).mkdir(parents=True, exist_ok=True)
        c = ROOT.TCanvas()
        stack.Draw("nostackb hist")
        self.apply_stack_properties(stack)
        self.apply_canvas_properties(c)
        if plot_dir:
            c.SaveAs(f"{plot_dir}/{plot_name}_no_legend.pdf")
        l.Draw()
        c.Draw()
        self.canvases[plot_name] = c
        if plot_dir:
            c.SaveAs(f"{plot_dir}/{plot_name}.pdf")