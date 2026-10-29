FEASIBILITY RECOVERY, NOT A CONVERGED OPTIMUM

This is the first converged SLSQP shape trial with a +0.01 degree incidence correction to move its drag below the fixed baseline cap. It is an independently evaluated feasible fallback design, not an optimizer termination result. Its run_complete.json intentionally records status -999 for independent evaluation only.

The top-level optimization is continuing separately. This folder preserves a feasible improvement if later trials fail. Frozen-connectivity results: CL=-1.1786084734, CD=0.0229375328, downforce increase 2.3265% relative to the same baseline, clearance 0.05018838 m, maximum thickness 0.05983021 m, y+ max 0.4421.

full_connectivity_check/ independently rebuilds the overset connectivity before CFD. No adjoint is evaluated there; the copied derivative-check file describes the parent frozen setup and does not validate unsupported full-mode derivatives. The full-connectivity solution also converges and is feasible: CL=-1.1782201460, CD=0.0229299363. This checks the frozen-connectivity approximation for this particular small deformation; it is not a grid-independence study.
