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


    def define_optimal_observables(self, names: str):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self.Define(f"O_{name}", f"{1/var} * (reco_sqme_12_nominal + reco_sqme_21_nominal - reco_sqme_12_{name} - reco_sqme_21_{name}) / (reco_sqme_12_nominal + reco_sqme_21_nominal)")
            self._define((f"mc_O_{name}", f"{1/var} * (mc_sqme_nominal - mc_sqme_{name}) / mc_sqme_nominal"), self._signal_categories)

