from analysis_framework.Analysis import Analysis
from AltSetupHandler import AltSetupHandler


class OptimalObservableHelper(Analysis):

    def __init__(self, dataset):
        super().__init__(dataset)
    

    # exists both here and in the reweight helper :/
    def define_optimal_observables_truth(self, names: list[str], truth_categories: list[str]):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self._define((f"mc_O_{name}", f"{1/var} * (mc_sqme_{name} - mc_sqme_nominal) / mc_sqme_nominal"), truth_categories)


    def define_optimal_observables_reco(self, names: list[str]):
        for name in names:
            var = AltSetupHandler.get_var_from_name_1d(name)
            self.Define(f"O_{name}", f"{1/var} * ((reco_sqme_12_{name} + reco_sqme_21_{name}) - (reco_sqme_12_nominal + reco_sqme_21_nominal)) / (reco_sqme_12_nominal + reco_sqme_21_nominal)")
