! Constant, family-specific physical wall velocities for steady ADflow.
! Velocities are prescribed data, not aerodynamic design variables.
module prescribedWallMotion
    use constants, only: intType, realType
    implicit none
    integer, parameter :: maxPrescribedWalls = 64
    integer(kind=intType), save :: wallFamilies(maxPrescribedWalls) = -1
    real(kind=realType), save :: wallVelocities(3, maxPrescribedWalls) = 0.0_realType
contains
    integer function wallVelocitySlot(family) result(slot)
        integer(kind=intType), intent(in) :: family
        integer :: k
        slot = 0
        do k = 1, maxPrescribedWalls
            if (wallFamilies(k) == family) then
                slot = k
                return
            end if
        end do
    end function wallVelocitySlot

    subroutine setWallVelocity(families, velocity, nFamilies, ierr)
        integer(kind=intType), intent(in) :: nFamilies, families(nFamilies)
        real(kind=realType), intent(in) :: velocity(3)
        integer(kind=intType), intent(out) :: ierr
        integer :: k, slot, freeSlot, j
        ierr = 0
        do k = 1, nFamilies
            slot = wallVelocitySlot(families(k))
            if (slot == 0) then
                freeSlot = 0
                do j = 1, maxPrescribedWalls
                    if (wallFamilies(j) < 0) then
                        freeSlot = j
                        exit
                    end if
                end do
                if (freeSlot == 0) then
                    ierr = 1
                    return
                end if
                slot = freeSlot
            end if
            wallFamilies(slot) = families(k)
            wallVelocities(:, slot) = velocity
        end do
    end subroutine setWallVelocity

    subroutine clearWallVelocities()
        wallFamilies = -1
        wallVelocities = 0.0_realType
    end subroutine clearWallVelocities

    subroutine getWallVelocityStats(families, nFamilies, sps, vmin, vmax, nFaces)
        use block, only: flowDoms, nDom
        use iteration, only: groundLevel
        use flowVarRefState, only: uRef
        integer(kind=intType), intent(in) :: nFamilies, families(nFamilies), sps
        real(kind=realType), intent(out) :: vmin(3), vmax(3)
        integer(kind=intType), intent(out) :: nFaces
        integer :: nn, mm, i, j, k
        real(kind=realType) :: velocity(3)
        vmin = huge(1.0_realType)
        vmax = -huge(1.0_realType)
        nFaces = 0
        do nn = 1, nDom
            do mm = 1, flowDoms(nn, groundLevel, sps)%nViscBocos
                if (.not. any(families == flowDoms(nn, groundLevel, sps)%BCData(mm)%famID)) cycle
                do j = flowDoms(nn, groundLevel, sps)%BCData(mm)%jcBeg, &
                       flowDoms(nn, groundLevel, sps)%BCData(mm)%jcEnd
                    do i = flowDoms(nn, groundLevel, sps)%BCData(mm)%icBeg, &
                           flowDoms(nn, groundLevel, sps)%BCData(mm)%icEnd
                        do k = 1, 3
                            velocity(k) = flowDoms(nn, groundLevel, sps)%BCData(mm)%uSlip(i, j, k) * uRef
                        end do
                        vmin = min(vmin, velocity)
                        vmax = max(vmax, velocity)
                        nFaces = nFaces + 1
                    end do
                end do
            end do
        end do
    end subroutine getWallVelocityStats
end module prescribedWallMotion
