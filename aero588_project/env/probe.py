import inspect, adflow, pyhyp, pygeo, idwarp, pyoptsparse, prefoil
from pygeo import DVGeometry, DVConstraints
for m in (adflow, pyhyp, pygeo, idwarp, pyoptsparse, prefoil):
    print(m.__name__, getattr(m, "__version__", "?"))
for n in ("addSpanwiseLocalDV", "addLocalDV", "addRefAxis", "addGlobalDV", "getLocalIndex"):
    print(n, hasattr(DVGeometry, n), str(inspect.signature(getattr(DVGeometry, n)))[:150] if hasattr(DVGeometry, n) else "")
for n in ("addLeTeConstraints", "addThicknessConstraints2D", "addVolumeConstraint", "addLinearConstraintsShape"):
    print(n, hasattr(DVConstraints, n), str(inspect.signature(getattr(DVConstraints, n)))[:150])
from pyoptsparse import OPT
for o in ("SLSQP", "IPOPT", "SNOPT"):
    try: OPT(o); print(o, "ok")
    except Exception as e: print(o, "no:", str(e)[:60])
