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


def build_sim_ws(runs, input_path, oo_name, run_ab, split_helicity_reversal: bool = True, pol_constraint: float | None = None):
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
        pars = fit_inputs[0]
        n_per_ab.append(fit_inputs[1])
        slopes.append(fit_inputs[3])
        Cov_O = fit_inputs[4]
        Cov_mean = np.asarray(Cov_O) / (n_per_ab[-1] * run_ab * lumi_share)
        nom_means = fit_inputs[2]
        if "cc10" in oo_name:
            # need to do some dirty removal of last par because the pol stuff in the cov matrix makes it singular
            pars = pars[:-1]
            Cov_mean = Cov_mean[:-1,:-1]
            nom_means = nom_means[:-1]
        Cov_mean_O.append(Cov_mean)
        pars_list.append(pars)
        nominal_means.append(nom_means)
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
        pars.remove("epol") if "epol" in pars else None
        pars.remove("ppol") if "ppol" in pars else None
        # pars = pars[:-2]  # remove pol pars
    fit_parameters = [ROOT.RooRealVar(par, par, 0., -0.5, 0.5) for par in pars]
    e_pol_fit_parameters = []
    p_pol_fit_parameters = []
    if split_helicity_reversal:
        e_pol_fit_parameters = [
            ROOT.RooRealVar("e_pol_L", "e_pol_L", 0., -0.5, 0.5),
            ROOT.RooRealVar("e_pol_R", "e_pol_R", 0., -0.5, 0.5),
            ROOT.RooRealVar("e_pol_0", "e_pol_0", 0., -0.5, 0.5)
            ]
        p_pol_fit_parameters = [
            ROOT.RooRealVar("p_pol_L", "p_pol_L", 0., -0.5, 0.5),
            ROOT.RooRealVar("p_pol_R", "p_pol_R", 0., -0.5, 0.5),
            ROOT.RooRealVar("p_pol_0", "p_pol_0", 0., -0.5, 0.5)
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
            if e_pol < 0:
                e_pol_par = e_pol_fit_parameters[0]  # L
            elif e_pol > 0:
                e_pol_par = e_pol_fit_parameters[1]  # R
            else:
                e_pol_par = e_pol_fit_parameters[2]  # 0
            if p_pol < 0:
                p_pol_par = p_pol_fit_parameters[0]  # L
            elif p_pol > 0:
                p_pol_par = p_pol_fit_parameters[1]  # R
            else:
                p_pol_par = p_pol_fit_parameters[2]  # 0
            run_fit_parameters += [e_pol_par, p_pol_par]
        else:
            e_pol_par = run_fit_parameters[-2]
            p_pol_par = run_fit_parameters[-1]
        build_model(w, run_fit_parameters, obs, nominal_means_run, slopes_run, Cov_mean_O_run, f"_run{i}")
        fit_model = w.pdf(f"fit_model_run{i}")
        if pol_constraint:
            # FIXME: this is wrong as the constraints are added multiple times :(
            # as a workaround we now figure out how often each constraint is applied
            # and scale up the value accordingly...
            # assume we only have the cases of 1 run or 4 runs (LL LR RL RR)
            n_same_constr = 1
            if len(runs) == 4:
                n_same_constr = 4
                if split_helicity_reversal:
                    n_same_constr /= 2
            scaled_pol_constraint = pol_constraint * np.sqrt(n_same_constr)
            e_pol_cnstr_v = abs(e_pol) * scaled_pol_constraint if e_pol != 0 else scaled_pol_constraint
            p_pol_cnstr_v = abs(p_pol) * scaled_pol_constraint if p_pol != 0 else scaled_pol_constraint
            e_pol_cnstr_pdf = ROOT.RooGaussian(f"constraint_e_pol_run{i}", f"constraint_e_pol_run{i}", e_pol_par, ROOT.RooFit.RooConst(0.), ROOT.RooFit.RooConst(e_pol_cnstr_v))
            p_pol_cnstr_pdf = ROOT.RooGaussian(f"constraint_p_pol_run{i}", f"constraint_p_pol_run{i}", p_pol_par, ROOT.RooFit.RooConst(0.), ROOT.RooFit.RooConst(p_pol_cnstr_v))
            w.Import(e_pol_cnstr_pdf)
            w.Import(p_pol_cnstr_pdf)
            fit_model = ROOT.RooProdPdf(f"fit_model_run{i}_with_constraints", f"fit_model_run{i}_with_constraints", [fit_model, w.pdf(e_pol_cnstr_pdf.GetName()), w.pdf(p_pol_cnstr_pdf.GetName())])
        models.append(fit_model)

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


def fit(workspace, model_name, ds, silent: bool = False, minos: bool = False, constant_parameters: list[str] = [], constraints: dict[str, float] = {}):
    """Fit the model to the given dataset and return the fitted parameters and their errors."""
    workspace.loadSnapshot("nominal_parameters")
    model = workspace.pdf(model_name)
    for par_name in constant_parameters:
        par = workspace.var(par_name)
        if par:
            par.setConstant(True)
        else:
            print(f"Warning: parameter {par_name} not found in workspace, cannot set it constant.")
    constraint_pdfs = []
    for par_name, sigma in constraints.items():
        par = workspace.var(par_name)
        if par:
            mean = par.getValV()
            gauss = ROOT.RooGaussian(f"constraint_{par_name}", f"constraint_{par_name}", par, ROOT.RooFit.RooConst(mean), ROOT.RooFit.RooConst(sigma))
            constraint_pdfs.append(gauss)
        else:
            print(f"Warning: parameter {par_name} not found in workspace, cannot add constraint.")
    if constraint_pdfs:
        all_pdfs = [model] + constraint_pdfs
        model = ROOT.RooProdPdf(f"{model_name}_with_constraints", f"{model_name}_with_constraints", all_pdfs)

    res = model.fitTo(ds, Minos=minos, Save=True, PrintLevel=-1 if silent else 1)
    return res


def make_par_histo(fit_res, parameter_names: list[str], x_label_dict: dict[str, str] = {}):
    n = len(parameter_names)
    h = ROOT.TH1D("", "", n, 0, n)
    fit_parameters = fit_res.floatParsFinal()
    for name in parameter_names:
        nice_name = x_label_dict[name] if name in x_label_dict else name
        try:
            par = fit_parameters.find(name)
            h.Fill(nice_name, par.getError())
        except:
            h.Fill(nice_name, 0)
    return h


def make_par_histos(fit_results, histos, stack, l, parameter_names: list[str],
                    legend_name_dict: dict[str, str] = {}, run_names: list[str] = [],
                    x_label_dict: dict[str, str] = {}):
    for i, name, in enumerate(run_names):
        fit_res = fit_results[name]
        h = make_par_histo(fit_res, parameter_names, x_label_dict=x_label_dict)
        h.SetFillColor(ROOT.kP10Blue + i)
        h.SetLineColor(ROOT.kP10Blue + i)
        leg_name = legend_name_dict[name] if name in legend_name_dict else name
        l.AddEntry(h, leg_name, "f")
        histos[name] = h
        stack.Add(h)


def make_plots_runs_per_oo_name(oo_name: str, fit_results: dict, run_names: list[str],
                                parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"],
                                legend_pars = None, legend_name_dict: dict[str, str] = {},
                                x_label_dict: dict[str, str] = {}):
    histos = {}
    stack = ROOT.THStack()
    if legend_pars:
        l = ROOT.TLegend(*legend_pars)
    else:
        l = ROOT.TLegend()
    l.SetHeader(legend_name_dict[oo_name] if oo_name in legend_name_dict else oo_name)
    make_par_histos(fit_results[oo_name], histos, stack, l, parameter_names, legend_name_dict, run_names, x_label_dict=x_label_dict)
    return histos, stack, l


def make_plots_oo_names_per_run_name(run_name: str, fit_results: dict, oo_names: list[str],
                                     parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"],
                                     legend_pars = None, legend_name_dict: dict[str, str] = {},
                                     x_label_dict: dict[str, str] = {}):
    histos = {}
    stack = ROOT.THStack()
    if legend_pars:
        l = ROOT.TLegend(*legend_pars)
    else:
        l = ROOT.TLegend()
    l.SetHeader(legend_name_dict[run_name] if run_name in legend_name_dict else run_name)
    for i, oo_name in enumerate(oo_names):
        fit_res = fit_results[oo_name][run_name]
        h = make_par_histo(fit_res, parameter_names, x_label_dict=x_label_dict)
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


    @staticmethod
    def apply_canvas_properties(c):
        c.SetLeftMargin(0.15)
        c.SetRightMargin(0.05)
        # c.SetBottomMargin(0.1375)
        c.SetBottomMargin(0.14)

    @staticmethod
    def apply_stack_properties(stack, unit:str = ""):
        # stack.GetXaxis().SetTitleSize(0.5)
        stack.GetXaxis().SetLabelSize(0.125)
        if unit:
            stack.GetYaxis().SetTitle(f"absolute uncertainty [{unit}]")
        else:
            stack.GetYaxis().SetTitle("absolute uncertainty")


    def draw_plots_runs_per_oo_name(self, plot_name: str, oo_name: str, fit_results: dict, run_names: list[str],
                                parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"],
                                legend_pars = None, legend_name_dict: dict[str, str] = {}, plot_dir: str|None = None, no_draw=False,
                                x_label_dict: dict[str, str] = {}, unit: str = ""):
        histos, stack, l = make_plots_runs_per_oo_name(oo_name, fit_results, run_names, parameter_names, legend_pars, legend_name_dict, x_label_dict=x_label_dict)
        self.histos[plot_name] = histos
        self.stacks[plot_name] = stack
        self.legends[plot_name] = l
        c = ROOT.TCanvas()
        stack.Draw("nostackb hist")
        self.apply_stack_properties(stack, unit=unit)
        self.apply_canvas_properties(c)
        c.Update()
        # stack.Draw("nostackb hist")
        if plot_dir:
            c.SaveAs(f"{plot_dir}/{plot_name}_no_legend.pdf")
        l.Draw()
        if not no_draw:
            c.Draw()
        self.canvases[plot_name] = c
        if plot_dir:
            c.SaveAs(f"{plot_dir}/{plot_name}.pdf")


    def draw_plots_oo_names_per_run_name(self, plot_name: str, run_name: str, fit_results: dict, oo_names: list[str],
                                parameter_names: list[str] = ["g1z", "ka", "la", "mW", "e_pol_L", "e_pol_R", "p_pol_L", "p_pol_R"],
                                legend_pars = None, legend_name_dict: dict[str, str] = {}, plot_dir: str|None = None, no_draw=False,
                                x_label_dict: dict[str, str] = {}, unit: str = ""):
        histos, stack, l = make_plots_oo_names_per_run_name(run_name, fit_results, oo_names, parameter_names, legend_pars, legend_name_dict, x_label_dict=x_label_dict)
        self.histos[plot_name] = histos
        self.stacks[plot_name] = stack
        self.legends[plot_name] = l
        if plot_dir:
            Path(plot_dir).mkdir(parents=True, exist_ok=True)
        c = ROOT.TCanvas()
        stack.Draw("nostackb hist")
        self.apply_stack_properties(stack, unit=unit)
        self.apply_canvas_properties(c)
        if plot_dir:
            c.SaveAs(f"{plot_dir}/{plot_name}_no_legend.pdf")
        l.Draw()
        if not no_draw:
            c.Draw()
        self.canvases[plot_name] = c
        if plot_dir:
            c.SaveAs(f"{plot_dir}/{plot_name}.pdf")


# not in plotter class
def draw_overlaid_stacks(plot_name: str, stack_name: str, p_names: list[str], plotters: dict[str, Plotter],
                         plot_dir: str|None = None, legend_pos: tuple[float, float, float, float]|None = None,
                         y_min: float|None = None, y_max: float|None = None, legend_columns: int|None = None,
                         extra_legend_texts: list[str] = []):
    c = ROOT.TCanvas()
    l = plotters[p_names[1]].legends[stack_name]
    if legend_pos:
        l.SetX1NDC(legend_pos[0])
        l.SetY1NDC(legend_pos[1])
        l.SetX2NDC(legend_pos[2])
        l.SetY2NDC(legend_pos[3])
    if legend_columns:
        l.SetNColumns(legend_columns)

    for i, name in enumerate(p_names):
        p = plotters[name]
        histos = p.histos[stack_name].values()
        for j, h in enumerate(histos):
            if i == 0:
                # h.SetFillColorAlpha(h.GetFillColor(), 0.75)
                h.SetFillColorAlpha(h.GetFillColor(), 0.5)
                h.SetLineColor(h.GetFillColor())
            elif i == 2:
                # h.SetFillColor(ROOT.kGray+3)
                h.SetFillColor(ROOT.kBlack)
                # h.SetFillColorAlpha(ROOT.kGray+3, 0.2)
                # h.SetFillColor(ROOT.kWhite)
                # h.SetFillStyle(3544)
                # h.SetFillStyle(3013)
                h.SetFillStyle(3006)
                # h.SetFillStyle(1001)
                # h.SetLineColor(ROOT.kBlack)
                # h.SetLineStyle(ROOT.kSolid)
            if j == 0 and extra_legend_texts:
                l.AddEntry(h, extra_legend_texts[i], "f")
        s = p.stacks[stack_name]
        if i == 0:
            if y_min is not None:
                s.SetMinimum(y_min)
            if y_max is not None:
                s.SetMaximum(y_max)
            s.Draw("nostackb hist")
        else:
            s.Draw("nostackb same hist")
    c.SetGridy()
    Plotter.apply_canvas_properties(c)
    if plot_dir:
        c.SaveAs(f"{plot_dir}/{plot_name}_overlaid_stacks_no_legend.pdf")
    if l:
        l.Draw()
    c.Draw()
    if plot_dir:
        c.SaveAs(f"{plot_dir}/{plot_name}_overlaid_stacks.pdf")
    return c

def make_corr_mat(fit_res, parameter_names: list[str], w, name_dict: dict[str, str] = {}):
    h = ROOT.TH2D()
    for  p_name1 in parameter_names:
        p1 = w.var(p_name1)
        p1_nice_name = name_dict[p_name1] if p_name1 in name_dict else p_name1
        for  p_name2 in reversed(parameter_names):
            p2 = w.var(p_name2)
            val = fit_res.correlation(p1, p2)
            p2_nice_name = name_dict[p_name2] if p_name2 in name_dict else p_name2
            h.Fill(p1_nice_name, p2_nice_name, round(val, 4))
    return h

# class precision_calculator():
#     # dirty hardcoded pars
#     pars = ["g1z", "ka", "la", "mW", "e_pol_0", "p_pol_0"]

#     def __init__(self, path: str, n: int):
#         self._get_workspace(path)
#         self.n = n
#         tmp_cov = self.ws.pdf("fit_model_run0").covarianceMatrix()
#         n_pars = tmp_cov.GetNrows()
#         self.full_cov = np.asarray([[tmp_cov(i,j) for j in range(n_pars)] for i in range(n_pars)])
#         self.full_cov *= n


#     def _get_workspace(self, path: str):
#         with ROOT.TFile.Open(path) as f:
#             self.ws = f.Get("w")


#     def get_exp_prec(self, pars: list[str]):
#         idcs = [self.pars.index(par) for par in pars]
#         red_cov = self.full_cov[np.ix_(idcs, idcs)]
#         return np.sqrt(np.diag(np.linalg.inv(red_cov)) / self.n)
