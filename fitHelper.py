import ROOT
from ast import literal_eval
from itertools import combinations_with_replacement
import numpy as np


def load_fit_input(input_path):
    with open(input_path) as infile:
        lines = infile.readlines()
    n_per_ab = literal_eval(lines[1])
    nominal_means = literal_eval(lines[2])
    slopes = literal_eval(lines[3])
    Cov_O = literal_eval(lines[4])
    means_list = literal_eval(lines[5])
    print(lines)
    return n_per_ab, nominal_means, slopes, Cov_O, means_list


def create_fit_model(model_name, workspace, coupling_names, nominal_means, slopes, Cov_mean_O):
    couplings = []
    for i in range(3):
        cpl = ROOT.RooRealVar(coupling_names[i], coupling_names[i], 0., -0.5, 0.5)
        couplings.append(cpl)

    obs = []
    mus_exp = []
    all_vars = []
    rel_changes = []
    for i, m in enumerate(nominal_means):
        sigma = np.sqrt(Cov_mean_O[i,i])
        o = ROOT.RooRealVar(f"O_{i+1}", f"O_{i+1}", m, m-5*sigma, m+5*sigma)
        obs.append(o)

        vars = []
        for j in range(3):
            # var = ROOT.RooProduct(f"var_{i}_{j}", f"var_{i}_{j}", couplings[j], ROOT.RooFit.RooConst(slopes[j][i]))
            var = ROOT.RooProduct(f"var_{i}_{j}", f"var_{i}_{j}", couplings[j], ROOT.RooFit.RooConst(slopes[i][j]))
            vars.append(var)
        all_vars += vars
        nominal_exp = ROOT.RooFit.RooConst(m)
        # mu_exp = ROOT.RooAddition(f"mu_exp_{i}", f"mu_exp_{i}", [nominal_exp] + vars)
        rel_change = ROOT.RooAddition(f"rel_change_{i}", f"rel_change_{i}", [ROOT.RooFit.RooConst(1.0)] + vars)
        rel_changes.append(rel_change)
        mu_exp = ROOT.RooProduct(f"mu_exp_{i}", f"mu_exp_{i}", rel_change, nominal_exp)
        mus_exp.append(mu_exp)

    # print(all_vars)

    n_obs = len(nominal_means)
    Cov_mean_O_root = ROOT.TMatrixDSym(n_obs)
    for i, j in combinations_with_replacement(range(3), r=2):
        Cov_mean_O_root[i][j] = Cov_mean_O[i, j]
        Cov_mean_O_root[j][i] = Cov_mean_O[i, j]

    model = ROOT.RooMultiVarGaussian(model_name, model_name, obs, mus_exp, Cov_mean_O_root)
    workspace.Import(model)
    workspace.saveSnapshot("nominal_parameters", couplings)
    workspace.saveSnapshot("nominal_observables", obs)
    workspace.defineSet("observables", obs)
    workspace.defineSet("parameters", couplings)


def make_datasets(workspace, names, means_list):
    obs = workspace.set("observables")
    for name in names:
        means = means_list[name]
        print(obs)
        for i in range(len(obs)):
            obs[f"O_{i+1}"] = means[i]
            # print(obs[f"O_{i+1}"])
        ds = ROOT.RooDataSet(name, name, obs)
        ds.add(obs)
        workspace.Import(ds)


def fit(workspace, model_name, ds_name):
    workspace.loadSnapshot("nominal_parameters")
    model = workspace.pdf(model_name)
    ds = workspace.data(ds_name)
    model.fitTo(ds)
    params = workspace.set("parameters")
    for p in params:
        print(p.GetName(), p.getVal(), p.getError())