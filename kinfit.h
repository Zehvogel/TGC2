#define MARLIN_USE_ROOT
#include <JetFitObject.h>
#include <MomentumConstraint.h>
#include <MassConstraint.h>
#include <SoftGaussMassConstraint.h>
#include <OPALFitterGSL.h>
#include <NewFitterGSL.h>
#include <Math/Vector4D.h>

struct enuWFit {
    // using lvec = ROOT::Math::PxPyPzMVector;
    using lvec = ROOT::Math::PxPyPzEVector;

    struct FitResult {
        double prob;
        double chi2;
        int error;
        int dof;
        int iterations;
        lvec obj1;
        lvec obj2;
        lvec obj3;
        lvec obj4;
    };

    double x_angle, e_cms, m_W, width_W;
    std::vector<double> errE, errTheta, errPhi;

    // FIXME: errors for the two actual jets (obj3, obj4) should be different then for the e (obj1) and nu (obj2)
    enuWFit(std::vector<double> eErr, std::vector<double> thetaErr, std::vector<double> phiErr, double xAngle, double eCMS, double mW, double width_W)
        : errE(eErr), errTheta(thetaErr), errPhi(phiErr), x_angle(xAngle), e_cms(eCMS), m_W(mW), width_W(width_W) {}

    FitResult operator()(lvec in_obj1, lvec in_obj2, lvec in_obj3, lvec in_obj4, bool doFitFlag) {
        if (!doFitFlag) {
            FitResult result;
            result.prob = -1.;
            result.chi2 = -1.;
            result.error = -1;
            result.dof = -1;
            result.iterations = -1;
            result.obj1 = in_obj1;
            result.obj2 = in_obj2;
            result.obj3 = in_obj3;
            result.obj4 = in_obj4;
            return result;
        }
        return doFit(in_obj1, in_obj2, in_obj3, in_obj4);
    }

    FitResult doFit(lvec in_obj1, lvec in_obj2, lvec in_obj3, lvec in_obj4) {
        JetFitObject f_obj1(in_obj1.P(), in_obj1.Theta(), in_obj1.Phi(), errE[0], errTheta[0], errPhi[0], 0.);
        JetFitObject f_obj2(in_obj2.P(), in_obj2.Theta(), in_obj2.Phi(), errE[1], errTheta[1], errPhi[1], 0.);
        JetFitObject f_obj3(in_obj3.P(), in_obj3.Theta(), in_obj3.Phi(), errE[2], errTheta[2], errPhi[2], 0.);
        JetFitObject f_obj4(in_obj4.P(), in_obj4.Theta(), in_obj4.Phi(), errE[3], errTheta[3], errPhi[3], 0.);

        MomentumConstraint e_constraint(1.0, 0.0, 0.0, 0.0, e_cms);
        MomentumConstraint px_constraint(0.0, 1.0, 0.0, 0.0, std::sin(x_angle/2.) * e_cms);
        MomentumConstraint py_constraint(0.0, 0.0, 1.0, 0.0, 0.0);
        MomentumConstraint pz_constraint(0.0, 0.0, 0.0, 1.0, 0.0);

        e_constraint.addToFOList(f_obj1);
        e_constraint.addToFOList(f_obj2);
        e_constraint.addToFOList(f_obj3);
        e_constraint.addToFOList(f_obj4);
        px_constraint.addToFOList(f_obj1);
        px_constraint.addToFOList(f_obj2);
        px_constraint.addToFOList(f_obj3);
        px_constraint.addToFOList(f_obj4);
        py_constraint.addToFOList(f_obj1);
        py_constraint.addToFOList(f_obj2);
        py_constraint.addToFOList(f_obj3);
        py_constraint.addToFOList(f_obj4);
        pz_constraint.addToFOList(f_obj1);
        pz_constraint.addToFOList(f_obj2);
        pz_constraint.addToFOList(f_obj3);
        pz_constraint.addToFOList(f_obj4);

        // MassConstraint w_constraint(m_W);
        SoftGaussMassConstraint w_constraint(width_W/2., m_W);

        w_constraint.addToFOList(f_obj3);
        w_constraint.addToFOList(f_obj4);

        // MassConstraint w_constraint2(m_W);
        SoftGaussMassConstraint w_constraint2(width_W/2., m_W);
        w_constraint2.addToFOList(f_obj1);
        w_constraint2.addToFOList(f_obj2);

        OPALFitterGSL fitter;
        // NewFitterGSL fitter;
        fitter.addFitObject(f_obj1);
        fitter.addFitObject(f_obj2);
        fitter.addFitObject(f_obj3);
        fitter.addFitObject(f_obj4);
        fitter.addConstraint(e_constraint);
        fitter.addConstraint(px_constraint);
        fitter.addConstraint(py_constraint);
        fitter.addConstraint(pz_constraint);
        fitter.addConstraint(w_constraint);
        // fitter.addConstraint(w_constraint2);

        // fitter.initialize();
        // fitter.setDebug(4);
        fitter.fit();

        FitResult result;
        result.prob = fitter.getProbability();
        result.chi2 = fitter.getChi2();
        result.error = fitter.getError();
        result.dof = fitter.getDoF();
        result.iterations = fitter.getIterations();
        // result.obj1 = lvec(f_obj1.getPx(), f_obj1.getPy(), f_obj1.getPz(), f_obj1.getMass());
        // result.obj2 = lvec(f_obj2.getPx(), f_obj2.getPy(), f_obj2.getPz(), f_obj2.getMass());
        // result.obj3 = lvec(f_obj3.getPx(), f_obj3.getPy(), f_obj3.getPz(), f_obj3.getMass());
        // result.obj4 = lvec(f_obj4.getPx(), f_obj4.getPy(), f_obj4.getPz(), f_obj4.getMass());
        result.obj1 = lvec(f_obj1.getPx(), f_obj1.getPy(), f_obj1.getPz(), f_obj1.getE());
        result.obj2 = lvec(f_obj2.getPx(), f_obj2.getPy(), f_obj2.getPz(), f_obj2.getE());
        result.obj3 = lvec(f_obj3.getPx(), f_obj3.getPy(), f_obj3.getPz(), f_obj3.getE());
        result.obj4 = lvec(f_obj4.getPx(), f_obj4.getPy(), f_obj4.getPz(), f_obj4.getE());

        return result;
    }
};

