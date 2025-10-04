from analysis_framework.Analysis import Analysis
import ROOT
from OO.whizard.model_parser import ModelParser
import subprocess
from AltSetupHandler import AltSetupHandler


class ReweightingHelper(Analysis):

    _omega_wrappers: dict
    def __init__(self, dataset):
        self._omega_wrappers = {}
        super().__init__(dataset)


    def initialise_omega_wrappers(self, configurations: dict[str,dict[str, float]]):
        whizard_prefix = subprocess.run(['whizard-config', '--prefix'], capture_output=True, encoding='ascii').stdout.strip()
        whizard_libs = f"{whizard_prefix}/lib/"
        # print(whizard_libs)
        ROOT.gSystem.AddDynamicPath(whizard_libs)
        ROOT.gSystem.Load("libwhizard.so")
        ROOT.gSystem.Load("libwhizard_main.so")
        ROOT.gSystem.Load("libomega.so")
        ROOT.gSystem.Load("OO/whizard/cc20_ac_inclusive/.libs/default_lib.so")
        ROOT.gInterpreter.Declare("#include \"OO/whizard/OmegaWrapper.h\"")

        model_parser = ModelParser("OO/whizard/SM_ac.mdl")
        # add derivation of lz and kz according to lep parametrisation
        model_parser.add_derived_parameter("lz", "la")
        model_parser.add_derived_parameter("kz", "1.0 - (ka - 1.0) * sw**2/cw**2 + (g1z - 1.0)")
        self._omega_wrappers["nominal"] = ROOT.OmegaWrapper(model_parser.get_parameters_list())

        for name, pars in configurations.items():
            model_parser.set_parameters(pars)
            self._omega_wrappers[name] = ROOT.OmegaWrapper(model_parser.get_parameters_list())


    def book_weights(self, categories: list[str], columns: list[str]):
        momenta = [f"{col}.energy(), {col}.Px(), {col}.Py(), {col}.Pz()" for col in columns]
        self.define_only_on(categories, "mc_ME_flv", "true_lep_charge > 0 ? 1 : 2")
        self.define_only_on(categories, "mc_ME_momenta", f"""
                    std::vector<double>({{
                            {','.join(momenta)}
                    }})
                    """)
        for name, omw in self._omega_wrappers.items():
            self.define_only_on(categories, f"mc_sqme_{name}", omw, ["mc_ME_momenta", "mc_ME_flv"])
            # divide by recalculated nominal as all the ILD values are broken...
            self.define_only_on(categories, f"weight_{name}", f"mc_sqme_{name} / mc_sqme_nominal")


    def define_optimal_observables_truth(self, names: list[str], truth_categories: list[str]):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self._define((f"mc_O_{name}", f"{1/var} * (mc_sqme_{name} - mc_sqme_nominal) / mc_sqme_nominal"), truth_categories)


    def define_optimal_observables_reco(self, names: list[str]):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self.Define(f"O_{name}", f"{1/var} * ((reco_sqme_12_{name} + reco_sqme_21_{name}) - (reco_sqme_12_nominal + reco_sqme_21_nominal)) / (reco_sqme_12_nominal + reco_sqme_21_nominal)")


    # TODO: some kind of prefix etc. would be good to add everywhere when the variations also come in
    def calc_reco_sqme(self, columns: list[str]):
        momenta = [f"{col}.energy(), {col}.Px(), {col}.Py(), {col}.Pz()" for col in columns]
        momenta2 = momenta.copy()
        # switch the jets
        momenta2[-2], momenta2[-1] = momenta2[-1], momenta2[-2]
        self.Define("reco_ME_flv", "iso_lep_charge > 0 ? 1 : 2")
        self.Define("reco_ME_momenta_12", f"""
            std::vector<double>({{
                    125., 0., 0., 125.,
                    125., 0., 0., -125.,
                    {','.join(momenta)}
            }})
        """)
        self.Define("reco_ME_momenta_21", f"""
            std::vector<double>({{
                    125., 0., 0., 125.,
                    125., 0., 0., -125.,
                    {','.join(momenta2)}
            }})
        """)
        for name, omw in self._omega_wrappers.items():
            self.Define(f"reco_sqme_12_{name}", omw, ["reco_ME_momenta_12", "reco_ME_flv"])
            self.Define(f"reco_sqme_21_{name}", omw, ["reco_ME_momenta_21", "reco_ME_flv"])
