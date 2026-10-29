The first frozen-connectivity attempt was interrupted during the second shape trial, which stalled in approximate Newton near residual 2. The preceding trial converged, CL=-1.1806621844 and CD=0.0229797932, but exceeded the drag cap and was not a feasible optimized result.

The continued run uses this optimization.hst for pyOptSparse hot-start replay. Its drag cap is kept unchanged through source/replay_settings.json. The production source changes NKSwitchTol from 1e-6 to 1e-4 and retries a stalled solve at 1e-2. It also returns explicit finite penalties for failed CFD evaluations because the SLSQP wrapper ignores fail flags. The objective, variable parameterization, bounds, geometry constraints, and SLSQP settings are unchanged.

Files here are an incomplete attempt, not a converged optimization. Production results and final termination status are in the parent run folder. Its wall time refers to the continued driver; this attempt's timings are in this log and history.
