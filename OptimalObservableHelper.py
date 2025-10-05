from analysis_framework.Analysis import Analysis
from AltSetupHandler import AltSetupHandler
from itertools import combinations_with_replacement
import ROOT


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

