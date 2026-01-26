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
        self._define(("true_lep_charge", "MCParticlesSkimmed.PDG[true_lep_idx] > 0. ? -1. : 1."), categories)
        self._define(("true_beam_e_lvec", make_lvec_M("MCParticlesSkimmed", "4")), categories)
        self._define(("true_beam_p_lvec", make_lvec_M("MCParticlesSkimmed", "5")), categories)
        self._define(("true_isr1_lvec", make_lvec_M("MCParticlesSkimmed", "6")), categories)
        self._define(("true_isr2_lvec", make_lvec_M("MCParticlesSkimmed", "7")), categories)

        self._define(("true_quark1_postCRC_idx", "return StableArgsort(MCParticlesSkimmed.generatorStatus == 2 && abs(MCParticlesSkimmed.PDG) < 6, [] (const auto a, const auto b) {return a > b;} )[2]"), categories)
        self._define(("true_quark2_postCRC_idx", "true_quark1_postCRC_idx+1"), categories)

        self._define(("true_quark1_postCRC_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark1_postCRC_idx")), categories)
        self._define(("true_quark2_postCRC_lvec", make_lvec_M("MCParticlesSkimmed", "true_quark2_postCRC_idx")), categories)

        self._define(("true_hadronic_W_postCRC_lvec", "true_quark1_postCRC_lvec + true_quark2_postCRC_lvec"), categories)

        # self.truth_defined = True
        self.truth_categories = categories


    def define_deltas(self, name, lvec1, lvec2, categories):
        self._define((f"{name}_delta_E", f"{lvec1}.energy() - {lvec2}.energy()"), categories)
        self._define((f"{name}_delta_E_rel_sqrt", f"({lvec1}.energy() - {lvec2}.energy()) / std::pow({lvec1}.energy(), 1.5)"), categories)
        self._define((f"{name}_delta_theta", f"{lvec1}.Theta() - {lvec2}.Theta()"), categories)
        self._define((f"{name}_delta_phi", f"ROOT::Math::VectorUtil::DeltaPhi({lvec1}, {lvec2})"), categories)
        self._define((f"{name}_delta_P", f"{lvec1}.P() - {lvec2}.P()"), categories)
        self._define((f"{name}_delta_P_rel_sqrt", f"({lvec1}.P() - {lvec2}.P()) / std::pow({lvec1}.P(), 1.5)"), categories)
        self._define((f"{name}_delta_Pt", f"{lvec1}.Pt() - {lvec2}.Pt()"), categories)
        self._define((f"{name}_delta_Pt2", f"{lvec1}.Perp2() - {lvec2}.Perp2()"), categories)
        self._define((f"{name}_delta_Px", f"{lvec1}.Px() - {lvec2}.Px()"), categories)
        self._define((f"{name}_delta_Py", f"{lvec1}.Py() - {lvec2}.Py()"), categories)
        self._define((f"{name}_delta_Pxy", f"{lvec1}.Px()*{lvec1}.Py() - {lvec2}.Px()*{lvec2}.Py()"), categories)
        self._define((f"{name}_delta_PxPy", f"({lvec1}.Px()+{lvec1}.Py()) - ({lvec2}.Px()+{lvec2}.Py())"), categories)
        self._define((f"{name}_delta_PxPy2", f"({lvec1}.Px()*{lvec1}.Px()+{lvec1}.Py()*{lvec1}.Py()) - ({lvec2}.Px()*{lvec2}.Px()+{lvec2}.Py()*{lvec2}.Py())"), categories)
        self._define((f"{name}_delta_Pz", f"{lvec1}.Pz() - {lvec2}.Pz()"), categories)
        self._define((f"{name}_delta_P_scaled2", f"{name}_delta_P / std::pow({lvec2}.P(), 1.5)"), categories)
        self._define((f"{name}_delta_P_scaled1", f"{name}_delta_P / std::pow({lvec1}.P(), 1.5)"), categories)