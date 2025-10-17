#include "edm4hep/ReconstructedParticleData.h"
// #include "RVec.hxx"
#include "podio/ObjectID.h"

using edm4hep::ReconstructedParticleData;
using PFOVec = ROOT::VecOps::RVec<ReconstructedParticleData>;
using IDVec = ROOT::VecOps::RVec<podio::ObjectID>;
using ROOT::VecOps::RVec;

PFOVec remove_constituents(const PFOVec& orig, const IDVec& orig_particles, const PFOVec& PFO_col, const ROOT::RVecI blacklist) {
    PFOVec res;
    res.reserve(orig.size());
    for (const auto& pfo : orig) {
        ReconstructedParticleData new_pfo(pfo);

        // subtract the blacklist
        for (unsigned int i = pfo.particles_begin; i < pfo.particles_end; i++) {
            auto id = orig_particles[i];
            if (blacklist[id.index]) {
                const auto& particle = PFO_col[id.index];
                new_pfo.energy -= particle.energy;
                new_pfo.momentum.x -= particle.momentum.x;
                new_pfo.momentum.y -= particle.momentum.y;
                new_pfo.momentum.z -= particle.momentum.z;
            }
        }

        res.push_back(new_pfo);
    }
    return res;
}

RVec<int> subset_to_mask(const PFOVec& orig_col, const IDVec& subset_col) {
    RVec<int> res(orig_col.size(), 0);

    for (const auto& id : subset_col) {
        res[id.index] = 1;
    }

    return res;
}

RVec<int> mcp_mask_to_pfo_mask(const RVec<int>& mcp_mask, const PFOVec& pfo_col, const IDVec& from, const IDVec& to, const RVec<float>& weights)
{
    RVec<int> pfo_mask(pfo_col.size(), 0);
    for (std::size_t i = 0; i < weights.size(); i++) {
        int j = to[i].index;
        if (mcp_mask[j]) {
            auto weight = weights[i];
            float c_weight = float( int(weight) / 10000 ) / 1000.f;
            if (c_weight > 0.5) {
                pfo_mask[from[i].index] = 1;
            }
        }
    }
    return pfo_mask;
}
// for each RecoMCTruth relation
//  take j = _to.index and check if mcp_mask[j] == 1
//  if yes:
//  take weight, take cluster weight from weight, check if > 0.5
//  take k = _from.index, and mark pfo_mask[k] == 1