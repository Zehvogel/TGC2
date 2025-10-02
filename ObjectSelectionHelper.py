from analysis_framework.Analysis import Analysis
import ROOT

def make_lvec_M(column: str, idx: str):
    return f"""
    ROOT::Math::PxPyPzMVector(
        {column}.momentum.x[{idx}],
        {column}.momentum.y[{idx}],
        {column}.momentum.z[{idx}],
        {column}.mass[{idx}]
    )
    """

def make_lvec_E(column: str, idx: str):
    return f"""
    ROOT::Math::PxPyPzEVector(
        {column}.momentum.x[{idx}],
        {column}.momentum.y[{idx}],
        {column}.momentum.z[{idx}],
        {column}.energy[{idx}]
    )
    """

class ObjectSelectionHelper(Analysis):

    truth_categories: list[str]


    def __init__(self, dataset):
        self.truth_categories = []
        super().__init__(dataset)


    def define_truth_objects(self, categories: list[str]):
        self._define(("true_quark1_idx", "ArgMax(MCParticlesSkimmed.generatorStatus == 2 && abs(MCParticlesSkimmed.PDG) < 6)"), categories)
        self._define(("true_quark2_idx", "true_quark1_idx + 1"), categories)
        self._define(("true_lep_idx", "true_quark1_idx + 2"), categories)
        self._define(("true_nu_idx", "true_quark1_idx + 3"), categories)
        self._define(("true_lep_lvec", make_lvec_M("MCParticlesSkimmed", "true_lep_idx")), categories)
        self._define(("true_nu_lvec", make_lvec_M("MCParticlesSkimmed", "true_nu_idx")), categories)
        self._define(("true_quark1_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark1_idx")), categories)
        self._define(("true_quark2_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark2_idx")), categories)
        self._define(("true_leptonic_W_lvec", "true_lep_lvec + true_nu_lvec"), categories)
        self._define(("true_hadronic_W_lvec", "true_quark1_lvec + true_quark2_lvec"), categories)
        self._define(("true_iso_lep_charge", "MCParticlesSkimmed.PDG[true_lep_idx] > 0. ? -1. : 1."), categories)

        # self.truth_defined = True
        self.truth_categories = categories

    def plot_resolution(self, column_name: str):
        """Plot (unscaled) resolution"""