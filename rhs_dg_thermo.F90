!***********************************************************************
!
!     Subroutine rhs_dg_thermo()
!
!     DG spatial calculations
!
!     Written by Colton J. Conroy (08-24-2018)
!
!***********************************************************************
subroutine rhs_dg_thermo(n)

use precisions
use global
implicit none

integer(int_p) :: left,right,iL,iR,check,kk,n,Nt
integer(int_p), parameter :: test = 2
real(real_p), parameter :: mu_c = 0.10d0, s_low = 2.0d0
real(real_p) :: nx, ny, tx, ty, l, speed, eig_L, eig_R,nt_ck, T_bound
real(real_p) :: T_speed, T_eig_L, T_eig_R, Tc, Tpt, vel, Tr, denX, denY
real(real_p) :: fhxL, fhyL, fuL, fvL, fhL,fuyL, fvyl, fhR, dhxL, dhyL
real(real_p) :: fuyR, fvxR, fhxR, fhyR, fuR, fvR, taux, tauy, dhxR, dhyR
real(real_p) :: fvyR, fuxR, fvxL, fuxL, den, mag_uv, sf
real(real_p) :: ftyR, ftxR, ftyL, ftxL, ftL, ftR, q_nR, q_tR
real(real_p) :: dvxL, dvxR, duyL, duyR, dvL, dvR, duL, duR, Texp, Tbl
real(real_p) :: s_layer_x, s_layer_y, b_layer_x, b_layer_y, mag_bl
real(real_p) :: jump_u, jump_v, jump_dudy, jump_dvdx, s_dudy, s_dvdx
real(real_p) :: fwxL, fwxR, fwyL,fwyR, Ck, alph_ck, s_mag
real(real_p) ::  jac, invjac, u_in(3), v_in(3), s_xn(3), s_yn(3)
real(real_p) :: d(3), huL(nlpts), huR(nlpts), hvL(nlpts), hvR(nlpts)
real(real_p) :: htL(nlpts), htR(nlpts), htj(napts), cosR(nlpts), fdhdy(nlpts)
real(real_p) :: etaL(nlpts), etaR(nlpts), cosL(nlpts), fdhdx(nlpts)
real(real_p) :: lambda_L(3), lambda_R(3), cos_n(3), cosj(napts)
real(real_p) :: fhu(nlpts), fhv(nlpts), fht(nlpts), fh(nlpts),fwy(nlpts)
real(real_p) :: uL(nlpts), vl(nlpts), hL(nlpts), hR(nlpts),fwx(nlpts)
real(real_p) :: uR(nlpts), vR(nlpts), dL(nlpts), dR(nlpts), fwj(napts)
real(real_p) :: cdotn(nlpts), qL_n(nlpts), qL_T(nlpts), blj(napts)
real(real_p) :: rL_n(nlpts), rL_T(nlpts), fdvdx(nlpts), fdudy(nlpts)
real(real_p) :: dvdxR(nlpts), dvdxL(nlpts), dudyR(nlpts), dudyL(nlpts)
real(real_p) :: wL(nlpts), wR(nlpts), wj(napts), dwdxj(napts), dwdyj(napts)
real(real_p) :: Bl(ndof,nlpts), Br(ndof,nlpts), rj(napts), dudyj(napts)
real(real_p) :: dj(napts), uj(napts), vj(napts), etaj(napts), dhdxj(napts)
real(real_p) :: huj(napts), hvj(napts), muj(napts), suj(napts),svj(napts)
real(real_p) :: fuxj(napts), fuyj(napts), fvxj(napts), fvyj(napts), dhdyj(napts)
real(real_p) :: ftxj(napts), ftyj(napts), sqj(napts), Qj(napts), dvdxj(napts)
real(real_p) :: fhxj(napts), fhyj(napts), hj(napts), sxj(napts), syj(napts)
real(real_p) :: Ax(ndof,napts), Ay(ndof,napts), Brf(ndof,nlpts)



Tr = sig*eff/(rho*shc*eHe)
Tc = 200.0d0
Ck = T_crust/T_int
!dvdx = 0.0d0; dudy = 0.0d0;
dhdx = 0.0d0; dhdy = 0.0d0;
dwdx = 0.0d0; dwdy = 0.0d0;
alph_ck = abs(alph - 1.0d0);
!.....Loop over edges

do j = 1, nedges

    nx = EDGEnormals(j,1)   ! normal vector @ edges
    ny = EDGEnormals(j,2)
    tx = -ny                ! tangential vector @ edges
    ty =  nx

!....Left variables

    left = EDGEelems(j,1)

    do i = 1, 3
        if (ELEMedges(left,i) == j) then
            iL = i
            goto 10
        endif
    enddo

10  huL  = matmul(PHIedge(:,:,iL),hu(:,left,irk))
    hvL  = matmul(PHIedge(:,:,iL),hv(:,left,irk))
    htL  = matmul(PHIedge(:,:,iL),ht(:,left,irk))
    etaL = matmul(PHIedge(:,:,iL),eta(:,left,irk))
    dvdxL= matmul(PHIedge(:,:,iL),dvdx(:,left,irk))
    dudyL= matmul(PHIedge(:,:,iL),dudy(:,left,irk))
    wL   = matmul(PHIedge(:,:,iL),w(:,left,irk))

    do i = 1, 3
        kk   = ELEMnodes(left,i)
        d(i)     = Z(kk)
        cos_n(i) = Zcos(kk)
    enddo

    dL   = matmul(PSIedge(:,:,iL),d)
    cosL = matmul(PSIedge(:,:,iL),cos_n)

    do i = 1, nlpts
        cdotn(i) = huL(i)*nx + hvL(i)*ny !...direction info for boundary flux
        qL_n(i)  = huL(i)*nx + hvL(i)*ny
        qL_T(i)  = huL(i)*tx + hvL(i)*ty
        rL_n(i)  = -huL(i)+2.0d0*sqrt(g/(dL(i)+etaL(i)))*etaL(i)
        rL_T(i)  = qL_T(i)
    enddo

!...Right variables and boundaries

    if (EDGEtype(j) == 0) then
        right = EDGEelems(j,2)
        do i = 1, 3
            if (ELEMedges(right,i) == j) then
                iR = i
                goto 20
            endif
        enddo

20      huR  = matmul(PHIedgef(:,:,iR),hu(:,right,irk))
        hvR  = matmul(PHIedgef(:,:,iR),hv(:,right,irk))
        htR  = matmul(PHIedgef(:,:,iR),ht(:,right,irk))
        etaR = matmul(PHIedgef(:,:,iR),eta(:,right,irk))
        dvdxR= matmul(PHIedgef(:,:,iR),dvdx(:,right,irk))
        dudyR= matmul(PHIedgef(:,:,iR),dudy(:,right,irk))
        wR   = matmul(PHIedgef(:,:,iR),w(:,right,irk))
        do i = 1, 3
            kk       = ELEMnodes(right,i)
            d(i)     = Z(kk)
            cos_n(i) = Zcos(kk)
        enddo

        dR    = matmul(PSIedge(:,:,iR),d)
        cosR  = matmul(PSIedge(:,:,iR),cos_n)

    elseif (EDGEtype(j) == 1) then             !...domain boundary

        if (inlet_edges(j) == 1) then !...inlet
            do i = 1, 3
                kk       = ELEMnodes(left,i)
                u_in(i)  = -u_inlet(kk)*1.0d0
                v_in(i)  = -v_inlet(kk)*1.0d0
            enddo
            uR  = matmul(PSIedge(:,:,iL),u_in)
            vR  = matmul(PSIedge(:,:,iL),v_in)
            hR  = dL + etaL
            do i = 1, nlpts
                huR(i)   = hR(i)*uR(i)
                hvR(i)   = hR(i)*vR(i)
                htR(i)   = hR(i)*T_int
                etaR(i)  = etaL(i)
                dR(i)    = dL(i)
                dvdxR(i) = dvdxL(i)
                dudyR(i) = dudyL(i)
                wR(i)    = wL(i)
            enddo
        elseif (outlet_edges(j) == 1) then !...outlet
            do i = 1, nlpts
                if (outlet_bc_type == 0) then
                    huR(i)  = huL(i)
                    hvR(i)  = hvL(i)
                elseif (outlet_bc_type == 1) then       !...radiation bc (NOT thermal radiation)
                    den    = 1.0d0/(nx*ty - ny*tx)
                    huR(i) = ( ty*rL_n(i) - ny*rL_T(i))*den
                    hvR(i) = (-tx*rL_n(i) + nx*rL_T(i))*den
                endif
                htR(i)   = htL(i)
                etaR(i)  = etaL(i)
                dR(i)    = dL(i)
                dvdxR(i) = dvdxL(i)
                dudyR(i) = dudyL(i)
                wR(i)    = wL(i)
            enddo
        else !...wall boundary
            do i = 1, nlpts
                    if (wall_bc_type == 0) then         !...no-slip
                        huR(i)  = 0.0d0
                        hvR(i)  = 0.0d0
                    elseif (wall_bc_type == 1) then     !...no-normal flow
                        den    = 1.0d0/(nx*ty - ny*tx)
                        q_nR   = -qL_n(i)
                        q_tR   =  qL_t(i)
                        huR(i) = ( ty*q_nR - ny*q_tR)*den
                        hvR(i) = (-tx*q_nR + nx*q_tR)*den
                    endif
                    htR(i)   = (dL(i) + etaL(i))*T_wall
                    etaR(i)  = etaL(i)
                    dR(i)    = dL(i)
                    dvdxR(i) = dvdxL(i)
                    dudyR(i) = dudyL(i)
                    wR(i)    = 0.0d0                     !...check this
            enddo
        endif
    endif

!...Numerical flux (Local Lax Friedrichs flux)

    do i = 1, nlpts

        hL(i) = dL(i) + etaL(i)
        hR(i) = dR(i) + etaR(i)

        if (hL(i) <= 1.0d-12) then
            uL(i)  = 0.0d0
            vL(i)  = 0.0d0
        else
            uL(i) = huL(i)/hL(i)
            vL(i) = hvL(i)/hL(i)
        endif

        if (hR(i) <= 1.0d-12) then
            uR(i) = 0.0d0
            vR(i) = 0.0d0
        else
            uR(i) = huR(i)/hR(i)
            vR(i) = hvR(i)/hR(i)
        endif
        fuxL   = huL(i)*uL(i) + half*g*(hL(i)**2.0d0 - dL(i)**2.0d0)*cosL(i)
        fuyL   = huL(i)*vL(i)
        fvxL   = hvL(i)*uL(i)
        fvyL   = hvL(i)*vL(i) + half*g*(hL(i)**2.0d0 - dL(i)**2.0d0)*cosL(i)
        ftxL   = htL(i)*uL(i)
        ftyL   = htL(i)*vL(i)
        fhxL   = hL(i)*uL(i)
        fhyL   = hL(i)*vL(i)
        dvxL   = vL(i)
        duyL   = uL(i)
        duL    = duyL*ny
        dvL    = dvxL*nx
        dhxL   = etaL(i)*nx
        dhyL   = etaL(i)*ny
        fwxL   = wL(i)*nx
        fwyL   = wL(i)*ny
        fuL    = fuxL*nx + fuyL*ny
        fvL    = fvxL*nx + fvyL*ny
        ftL    = ftxL*nx + ftyL*ny
        fhL    = fhxL*nx + fhyL*ny

        fuxR   = huR(i)*uR(i) + half*g*(hR(i)**2.0d0 - dR(i)**2.0d0)*cosR(i)
        fuyR   = huR(i)*vR(i)
        fvxR   = hvR(i)*uR(i)
        fvyR   = hvR(i)*vR(i) + half*g*(hR(i)**2.0d0 - dR(i)**2.0d0)*cosR(i)
        ftxR   = htR(i)*uR(i)
        ftyR   = htR(i)*vR(i)
        fhxR   = hR(i)*uR(i)
        fhyR   = hR(i)*vR(i)
        dvxR   = vR(i)
        duyR   = uR(i)
        dvR    = dvxR*nx
        duR    = duyR*ny
        dhxR   = etaR(i)*nx
        dhyR   = etaR(i)*ny
        fwxR   = wR(i)*nx
        fwyR   = wR(i)*ny
        fuR    = fuxR*nx + fuyR*ny
        fvR    = fvxR*nx + fvyR*ny
        ftR    = ftxR*nx + ftyR*ny
        fhR    = fhxR*nx + fhyR*ny


        lambda_L(1) = abs((uL(i)-sqrt(g*hL(i)*cosL(i)))*nx + (vL(i)-sqrt(g*hL(i)*cosL(i)))*ny)
        lambda_L(2) = abs(uL(i)*nx + vL(i)*ny)
        lambda_L(3) = abs((uL(i)+sqrt(g*hL(i)*cosL(i)))*nx + (vL(i)+sqrt(g*hL(i)*cosL(i)))*ny)
        lambda_R(1) = abs((uR(i)-sqrt(g*hR(i)*cosR(i)))*nx + (vR(i)-sqrt(g*hR(i)*cosR(i)))*ny)
        lambda_R(2) = abs(uR(i)*nx + vR(i)*ny)
        lambda_R(3) = abs((uR(i)+sqrt(g*hR(i)*cosR(i)))*nx + (vR(i)+sqrt(g*hR(i)*cosR(i)))*ny)

        eig_L = maxval(lambda_L);
        eig_R = maxval(lambda_R);
        speed = max(eig_L,eig_R);

        fhu(i) = half*((fuL+fuR)-abs(speed)*(huR(i)-huL(i)))
        fhv(i) = half*((fvL+fvR)-abs(speed)*(hvR(i)-hvL(i)))
        fht(i) = half*((ftL+ftR)-abs(speed)*(htR(i)-htL(i)))
        fh(i)  = half*((fhL+fhR)-abs(speed)*(hR(i)-hL(i)))

        fdvdx(i) = half*(dvL+dvR)
        fdudy(i) = half*(duL+duR)
        fdhdx(i) = half*(dhxL + dhxR)
        fdhdy(i) = half*(dhyL + dhyR)
        fwx(i)   = half*(fwxL+fwxR)
        fwy(i)   = half*(fwyL+fwyR)

    enddo

    Bl(:,:) = half*EDGElengths(j)*B(:,:,iL)
    rhs_hu(:,left,irk)  = rhs_hu(:,left,irk)  - matmul(Bl,fhu)
    rhs_hv(:,left,irk)  = rhs_hv(:,left,irk)  - matmul(Bl,fhv)
    rhs_ht(:,left,irk)  = rhs_ht(:,left,irk)  - matmul(Bl,fht)
    rhs_eta(:,left,irk) = rhs_eta(:,left,irk) - matmul(Bl,fh)
    dvdx(:,left,irk)    = dvdx(:,left,irk)    - matmul(Bl,fdvdx)
    dudy(:,left,irk)    = dudy(:,left,irk)    - matmul(Bl,fdudy)
    dhdx(:,left)        = dhdx(:,left)        - matmul(Bl,fdhdx)
    dhdy(:,left)        = dhdy(:,left)        - matmul(Bl,fdhdy)
    dwdx(:,left)        = dwdx(:,left)        - matmul(Bl,fwx)
    dwdy(:,left)        = dwdy(:,left)        - matmul(Bl,fwy)
    if (EDGEtype(j) == 0) then  !...Interior edges
        Br(:,:) = half*EDGElengths(j)*B(:,:,iR)
        do i = 1, nlpts
            do kk = 1, ndof
                Brf(kk,i) = Br(kk,nlpts+1-i)
            enddo
        enddo
        rhs_hu(:,right,irk)  = rhs_hu(:,right,irk)  + matmul(Brf,fhu)
        rhs_hv(:,right,irk)  = rhs_hv(:,right,irk)  + matmul(Brf,fhv)
        rhs_ht(:,right,irk)  = rhs_ht(:,right,irk)  + matmul(Brf,fht)
        rhs_eta(:,right,irk) = rhs_eta(:,right,irk) + matmul(Brf,fh)
        dvdx(:,right,irk)    = dvdx(:,right,irk)    + matmul(Brf,fdvdx)
        dudy(:,right,irk)    = dudy(:,right,irk)    + matmul(Brf,fdudy)
        dhdx(:,right)        = dhdx(:,right)        + matmul(Brf,fdhdx)
        dhdy(:,right)        = dhdy(:,right)        + matmul(Brf,fdhdy)
        dwdx(:,right)        = dwdx(:,right)        + matmul(Brf,fwx)
        dwdy(:,right)        = dwdy(:,right)        + matmul(Brf,fwy)
    endif
enddo


!...Element calculations

do j = 1, nelems

    do i = 1, 3
        kk = ELEMnodes(j,i)
        s_xn(i)  = Zx(kk)
        s_yn(i)  = Zy(kk)
        cos_n(i) = Zcos(kk)
        d(i)     = Z(kk)
    enddo

    sxj  = matmul(PSI,s_xn); sxj = g*sxj;
    syj  = matmul(PSI,s_yn); syj = g*syj;
    cosj = matmul(PSI,cos_n);
    dj   = matmul(PSI,d)

    huj  = matmul(PHI,hu(:,j,irk))
    hvj  = matmul(PHI,hv(:,j,irk))
    htj  = matmul(PHI,ht(:,j,irk))
    etaj = matmul(PHI,eta(:,j,irk))
    wj   = matmul(PHI,w(:,j,irk))

!...element areas and Jacobians
    Ax = half*(ELEMxy(j,4)*A(:,:,1)+ELEMxy(j,1)*A(:,:,2))
    Ay = half*(ELEMxy(j,3)*A(:,:,1)+ELEMxy(j,2)*A(:,:,2))
    jac    = ELEMarea(j)/2.0d0
    invjac = 2.0d0/ELEMarea(j)
!...calc gradients of H and w
    dhdx(:,j) = dhdx(:,j) + matmul(Ax,etaj)
    dhdy(:,j) = dhdy(:,j) + matmul(Ay,etaj)
    dhdx(:,j) = invjac*dhdx(:,j)
    dhdy(:,j) = invjac*dhdy(:,j)
    dwdx(:,j) = dwdx(:,j) + matmul(Ax,wj)
    dwdy(:,j) = dwdy(:,j) + matmul(Ay,wj)
    dwdx(:,j) = invjac*dwdx(:,j)
    dwdy(:,j) = invjac*dwdy(:,j)
    dwdxj= matmul(PHI,dwdx(:,j))
    dwdyj= matmul(PHI,dwdy(:,j))
    dhdxj= matmul(PHI,dhdx(:,j))
    dhdyj= matmul(PHI,dhdy(:,j))
    do i = 1, napts

        if (isnan(huj(i))) then
            write(*,*)'***** WARNING *****'
            write(*,*)'Solution is NaN @ time =',t,'secs.'
            write(*,*)'Execution will stop ...'
            write(*,*)'eta=',etaj
            call global_output()
            STOP
        endif

        s_mag = sqrt(sxj(i)**2.0d0 + syj(i)**2.0d0)

        hj(i) = dj(i) + etaj(i)

        if (hj(i) <= 1.0d-12) then
            uj(i)  = 0.0d0
            vj(i)  = 0.0d0
            Tpt    = T_int
            Tbl   = (Tpt - T_crust)*exp(-kt) + T_crust
            taux   = 0.0d0
            tauy   = 0.0d0
            sxj    = 0.0d0
            syj    = 0.0d0
            Qj(i)  = 0.0d0
        else
            uj(i) = huj(i)/hj(i)
            vj(i) = hvj(i)/hj(i)
            Tpt   = htj(i)/hj(i)
            Tbl   = (Tpt - T_crust)*exp(-kt) + T_crust
            vel   = sqrt(uj(i)**2.0d0 + vj(i)**2.0d0)

            if (vel > 1.0d-6 .and. Tpt > 1173.0d0) then        !...radiation
                Qj(i) = Tr/hj(i)*(Tpt**4.0d0 - T_air**4.0d0)
            elseif (vel > 1.0d-6 .and. Tpt <= 1173.0d0) then   !...convection
                Qj(i) = 0.50d0*Tc/hj(i)*(Tpt - T_air)**(4.0d0/3.0d0)
            else
                Qj(i) = 0.0d0
            endif
            mag_uv = sqrt(uj(i)**2.0d0 + vj(i)**2.0d0)
            denX = -dwdxj(i)
            if (abs(denX) < 1.0d-16) then
                s_layer_x = hj(i)*s_low
            else
                s_layer_x = abs(uj(i)/denX)
            endif
            if (abs(s_layer_x) < 1.0d-16) then
                s_layer_x = hj(i)*s_low
            endif
            denY = -dwdyj(i)
            if (abs(denY) < 1.0d-16) then
                s_layer_y = hj(i)*s_low
            else
                s_layer_y = abs(vj(i)/denY)
            endif
            if (abs(s_layer_y) < 1.0d-16) then
                s_layer_y = hj(i)*s_low
            endif
            b_layer_x = b_layer*s_layer_x
            b_layer_y = b_layer*s_layer_y
            blj(i)    = sqrt(b_layer_x**2.0d0 + b_layer_y**2.0d0)
            taux  = (uj(i)/(uj(i)*hj(i)*b_layer_x))*abs(uj(i)/(b_layer_x))**(alph-1.0d0)
            tauy  = (vj(i)/(vj(i)*hj(i)*b_layer_y))*abs(vj(i)/(b_layer_y))**(alph-1.0d0)
        endif
        if (Tbl <= 600.0d0) then
            write(*,*)'Temperature is too low'
            write(*,*)'T=',Tbl
            write(*,*)'Tpt=', Tpt
            write(*,*)'kt =', kt
            write(*,*)'Ck = ', Ck
            call global_output()
            stop
        endif
        Texp   = Amu + Bmu/(Tbl-Cmu)
        muj(i) = 10.0d0**Texp;
        muj(i) = muj(i)/(rho);
        taux   = muj(i)*taux
        tauy   = muj(i)*tauy
        if (alph_ck < 1.0d-12) then
            Texp   = Amu + Bmu/(Tpt-Cmu)
            muj(i) = 10.0d0**Texp;
            muj(i) = muj(i)/(rho);
        else
            Texp   = Amu + Bmu/(Tpt-Cmu)
            muj(i) = 10.0d0**Texp;
            mag_bl = sqrt(b_layer_x**2.0d0+b_layer_y**2.0d0)
            muj(i) = muj(i)*(mag_uv/mag_bl)**(alph-1.0d0)
            muj(i) = muj(i)/(rho);
        endif
        suj(i)    = -sxj(i)*hj(i) - taux*hj(i)*uj(i) - hj(i)*sign(tau_yield,uj(i))
        svj(i)    = -syj(i)*hj(i) - tauy*hj(i)*vj(i) - hj(i)*sign(tau_yield,vj(i))
        fuxj(i)   = huj(i)*uj(i) + half*g*(hj(i)**2.0d0 - dj(i)**2.0d0)*cosj(i)
        fuyj(i)   = huj(i)*vj(i)
        fvxj(i)   = hvj(i)*uj(i)
        fvyj(i)   = hvj(i)*vj(i) + half*g*(hj(i)**2.0d0 - dj(i)**2.0d0)*cosj(i)
        ftxj(i)   = htj(i)*uj(i)
        ftyj(i)   = htj(i)*vj(i)
        fhxj(i)   = huj(i)
        fhyj(i)   = hvj(i)
        fwj(i)    = -uj(i)*dhdxj(i) - vj(i)*dhdyj(i)
        dvdxj(i)  = vj(i)
        dudyj(i)  = uj(i)
    enddo

    mu(:,j)          = jac*matmul(C,muj)*rho; mu(:,j)  = invjac*mu(:,j);

    b_thick(:,j)     = jac*matmul(C,blj); b_thick(:,j) = invjac*b_thick(:,j);

    w(:,j,irk)       = jac*matmul(C,fwj); w(:,j,irk)   = invjac*w(:,j,irk);

    rhs_hu(:,j,irk)  = rhs_hu(:,j,irk)  + matmul(Ax,fuxj)   &
                                        + matmul(Ay,fuyj)   &
                                        + jac*matmul(C,suj)

    rhs_hv(:,j,irk)  = rhs_hv(:,j,irk)  + matmul(Ax,fvxj)   &
                                        + matmul(Ay,fvyj)   &
                                        + jac*matmul(C,svj)

    rhs_ht(:,j,irk)  = rhs_ht(:,j,irk)  + matmul(Ax,ftxj)   &
                                        + matmul(Ay,ftyj)   &
                                        + jac*matmul(C,Qj)

    rhs_eta(:,j,irk) = rhs_eta(:,j,irk) + matmul(Ax,fhxj)   &
                                        + matmul(Ay,fhyj)

    dvdx(:,j,irk)    = dvdx(:,j,irk)    + matmul(Ax,dvdxj)
    dudy(:,j,irk)    = dudy(:,j,irk)    + matmul(Ay,dudyj)

    rhs_hu(:,j,irk)  = invjac*rhs_hu(:,j,irk)
    rhs_hv(:,j,irk)  = invjac*rhs_hv(:,j,irk)
    rhs_ht(:,j,irk)  = invjac*rhs_ht(:,j,irk)
    rhs_eta(:,j,irk) = invjac*rhs_eta(:,j,irk)
    dvdx(:,j,irk)    = invjac*dvdx(:,j,irk)
    dudy(:,j,irk)    = invjac*dudy(:,j,irk)

enddo

end subroutine rhs_dg_thermo

