"""Apply a version-checked moving-wall extension to ADflow 2.12.2.

Generated forward/reverse sources are patched explicitly; rerunning Tapenade
would require regenerating these additions. Only constant wall data is supported.
"""
from pathlib import Path
import argparse
import difflib
import re
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--adflow-root', type=Path, default=Path('/home/mdolabuser/repos/adflow'))
args = parser.parse_args()
root = args.adflow_root
here = Path(__file__).resolve().parent
assert '__version__ = "2.12.2"' in (root/'adflow/__init__.py').read_text()
changes = []

def save(relative, old, new):
    assert old != new, relative
    p = root/relative
    p.write_text(new)
    changes.extend(difflib.unified_diff(old.splitlines(True), new.splitlines(True),
                                      fromfile='a/'+relative, tofile='b/'+relative))

module = (here/'prescribedWallMotion.F90').read_text()
save('src/modules/prescribedWallMotion.F90', '', module)
p = root/'src/build/fileList'
old = p.read_text()
anchor = '\t\tmodules/constants.F90\\\n'
assert old.count(anchor) == 1
save('src/build/fileList', old, old.replace(anchor, anchor+'\t\tmodules/prescribedWallMotion.F90\\\n'))

interface = '''
       module prescribedwallmotion
         use precision
         subroutine setwallvelocity(families,velocity,nfamilies,ierr)
           integer(kind=inttype), dimension(nfamilies),intent(in) :: families
           real(kind=realtype), dimension(3),intent(in) :: velocity
           integer(kind=inttype), optional,intent(in),check(len(families)>=nfamilies),depend(families) :: nfamilies=len(families)
           integer(kind=inttype),intent(out) :: ierr
         end subroutine setwallvelocity
         subroutine clearwallvelocities
         end subroutine clearwallvelocities
         subroutine getwallvelocitystats(families,nfamilies,sps,vmin,vmax,nfaces)
           integer(kind=inttype), dimension(nfamilies),intent(in) :: families
           integer(kind=inttype), optional,intent(in),check(len(families)>=nfamilies),depend(families) :: nfamilies=len(families)
           integer(kind=inttype),intent(in) :: sps
           real(kind=realtype), dimension(3),intent(out) :: vmin,vmax
           integer(kind=inttype),intent(out) :: nfaces
         end subroutine getwallvelocitystats
       end module prescribedwallmotion

'''
p = root/'src/f2py/adflow.pyf'
old = p.read_text()
anchor = '       module preprocessingapi\n'
assert old.count(anchor) == 1
save('src/f2py/adflow.pyf', old, old.replace(anchor, interface+anchor))

def override(derivative=False, seed_only=False):
    lines = ['    ! Prescribed absolute wall speed; constant with respect to shape/state.',
             '    block', '      integer :: mw, slot', '      do mw = 1, nViscBocos',
             '        slot = wallVelocitySlot(BCData(mw)%famID)', '        if (slot > 0) then']
    if not seed_only:
        for k in range(1, 4):
            lines.append(f'          BCData(mw)%uSlip(:, :, {k}) = wallVelocities({k}, slot) / uRef')
    if derivative or seed_only:
        lines.append('          BCDataD(mw)%uSlip = 0.0_realType')
    lines += ['        end if', '      end do', '    end block\n']
    return '\n'.join(lines)

sources = ['src/solver/solverUtils.F90',
           'src/adjoint/outputForward/solverUtils_d.f90',
           'src/adjoint/outputReverse/solverUtils_b.f90',
           'src/adjoint/outputReverseFast/solverUtils_fast_b.f90']
for relative in sources:
    old = (root/relative).read_text()
    pattern = re.compile(r'(?im)^([ \t]*)subroutine (slipvelocitiesfinelevel_block(?:_d|_b)?)\([^\n]*\).*?^\1end subroutine \2\b', re.S)
    seen = []
    def patch_function(match):
        name = match[2].lower()
        seen.append(name)
        code = match[0]
        code, n = re.subn(r'(?im)^(\s*use constants\s*)$',
                         r'\1\n    use prescribedWallMotion, only: wallVelocitySlot, wallVelocities', code, count=1)
        assert n == 1, name
        # The reverse seed of a prescribed constant must be discarded before
        # the existing reverse sweep propagates wall-motion derivatives.
        if name.endswith('_b'):
            marker = re.search(r'(?im)^\s*! determine the situation we are having here\.', code)
            assert marker, name
            code = code[:marker.start()] + '\n'+override(seed_only=True) + code[marker.start():]
        tail = re.search(r'(?im)^\s*end subroutine '+name+r'\b', code)
        code = code[:tail.start()]+'\n'+override(derivative=name.endswith('_d'))+code[tail.start():]
        return code
    new = pattern.sub(patch_function, old)
    expected = 1 if relative.endswith('solverUtils.F90') or 'ReverseFast' in relative else 2
    assert len(seen) == expected, (relative, seen)
    save(relative, old, new)

api = '''    def setWallVelocity(self, velocity, groupName):
        """Prescribe a constant physical wall velocity [m/s] for a viscous family.

        Steady, real-valued problems only. The caller must prescribe a tangential
        velocity on a fixed wall. Velocity is constant data, not a design variable.
        Call on every MPI rank before the flow solve. Airfoil families are untouched.
        """
        if self.getOption("equationMode").lower() != "steady":
            raise Error("Prescribed wall translation currently supports steady mode only")
        if numpy.issubdtype(numpy.dtype(self.dtype), numpy.complexfloating):
            raise Error("Prescribed wall translation requires the patched real solver")
        velocity = numpy.asarray(velocity, dtype=float)
        if velocity.shape != (3,) or not numpy.isfinite(velocity).all():
            raise Error("Wall velocity must be a finite three-component vector in m/s")
        families = numpy.asarray(self._getFamilyList(groupName), dtype=numpy.int32)
        if len(families) == 0:
            raise Error("The requested wall family group is empty")
        signature = (tuple(families.tolist()), tuple(velocity.tolist()))
        if any(item != signature for item in self.comm.allgather(signature)):
            raise Error("Every MPI rank must prescribe identical wall velocity data")
        ierr = self.adflow.prescribedwallmotion.setwallvelocity(families, velocity)
        if ierr:
            raise Error("Too many prescribed moving-wall families (maximum 64)")

    def getWallVelocityStats(self, groupName):
        """Return MPI-global min/max physical wall speeds and number of faces.

        This reads the actual internal no-slip boundary data after a flow solve.
        It does not return first-off-wall fluid velocities.
        """
        families = numpy.asarray(self._getFamilyList(groupName), dtype=numpy.int32)
        lo, hi, count = self.adflow.prescribedwallmotion.getwallvelocitystats(families, 1)
        count = self.comm.allreduce(int(count), op=MPI.SUM)
        if count == 0:
            raise Error("No viscous faces found in the requested wall family")
        globalLo, globalHi = numpy.empty_like(lo), numpy.empty_like(hi)
        self.comm.Allreduce(lo, globalLo, op=MPI.MIN)
        self.comm.Allreduce(hi, globalHi, op=MPI.MAX)
        return {"min_mps": globalLo.tolist(), "max_mps": globalHi.tolist(), "faces": count}

'''
p = root/'adflow/pyADflow.py'
old = p.read_text()
anchor = '    def setRotationRate('
assert old.count(anchor) == 1
save('adflow/pyADflow.py', old, old.replace(anchor, api+anchor))
(root/'adflow-moving-wall.patch').write_text(''.join(changes))
print('Patched ADflow 2.12.2: prescribed wall motion, Python/F2PY interfaces, and derivative routines.')
