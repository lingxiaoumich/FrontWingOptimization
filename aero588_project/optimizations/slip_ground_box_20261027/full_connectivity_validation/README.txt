INDEPENDENT FINAL CFD VALIDATION

This folder checks the selected final geometry with ADflow oversetUpdateMode=full, rebuilding hole cutting and interpolation. It starts from the selected design's converged flow states. No optimizer or adjoint is run here.

The independent solve converged and satisfies the surface box, minimum allowable maximum thickness, and original drag cap. The change in downforce relative to the frozen-connectivity optimization solution is -0.11543%. This comparison checks the interpolation treatment for the selected deformation; it is not a grid-independence study.

See validation_complete.json, final/result.json, settings.txt, and run.log. The parent comparison_results.txt includes both force predictions. The main before/after plots use the primary optimization-model solution.
