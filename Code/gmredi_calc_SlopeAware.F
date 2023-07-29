#include "GMREDI_OPTIONS.h"

CBOP
C !ROUTINE: SlopeAware

C !INTERFACE: ==========================================================
      SUBROUTINE GMREDI_CALC_SlopeAware(
     I             sigmaX, sigmaY, sigmaR,
     I             bi, bj, myTime, myIter, myThid, 
     O             SlopeAwareModi                    )

C     !DESCRIPTION:
C     *==========================================================*
C     | SUBROUTINE GMREDI_CALC_GEOM
C     | Calculate GM coefficient with the GEOMETRIC prescription
C     | GEOMK is located at the cell centre.
C     *==========================================================*

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "GMREDI.h"

C !INPUT PARAMETERS: ===================================================
C     sigmaXYR  :: density variables
C     bi, bj    :: tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number in simulation
C     myThid    :: My Thread Id. number

      INTEGER i,j,k
    
      _RL sigmaX(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER kp1
      
      _RL, save :: depthC(1-OLx:sNx+OLx, 1-OLy:sNy+OLy)
            
      _RL dHdy(1-Olx:sNx+Olx,1-Oly:sNy+Oly)  
      _RL SBurger(1-Olx:sNx+Olx,1-Oly:sNy+Oly) 
      _RL SDelta(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL SlopeAwareModi(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL N(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL N2(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL M2(1-Olx:sNx+Olx,1-Oly:sNy+Oly,Nr)
      _RL BVint(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL BVzavg(1-Olx:sNx+Olx,1-Oly:sNy+Oly)  
      _RL dSigmaDxZint(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL dSigmaDyZint(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL dSigmaDrZint(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL Siso(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL Beta_t(1-Olx:sNx+Olx,1-Oly:sNy+Oly)
      _RL Lrhines(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      
      
C       depths for doing depth-averaging
        depthC = 0.0 _d 0
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
              depthC(i,j) = depthC(i,j) + drF(k) * hFacC(i,j,k,bi,bj)
            ENDDO
          ENDDO
        ENDDO
        
      
      IF(GEOM_SloAwa_Burger) THEN

       DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
         BVint(i,j)  = 0. _d 0
        ENDDO
       ENDDO

C    Bathymetric slope at cell center
C    depthC used here is positive  
     
      DO j=1-Oly+1,sNy+Oly-1
       DO i=1-Olx,sNx+Olx
         dHdy(i,j) = -op5 *
     &         (maskS(i,j,1,bi,bj)
     &         *_recip_dyC(i,j,bi,bj)*
     &           (depthC(i,j)-depthC(i,j-1))
     &          +maskS(i,j+1,1,bi,bj)
     &         *_recip_dyC(i,j+1,bi,bj)*
     &           (depthC(i,j+1)-depthC(i,j)) )
       ENDDO
      ENDDO  
    
      
      


C     Square of the buoyancy frequency at the top of a grid cell
C     Enforce a minimum N2; Mask N2, so it is zero at bottom
C     Calculation of N2 and N follows gmredi_k3d.F
      DO k=2,Nr
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
         N2(i,j,k) = -gravity*recip_rhoConst*sigmaR(i,j,k)  
         N2(i,j,k) = MAX(N2(i,j,k),1.0 _d -8)*maskC(i,j,k,bi,bj)
         N(i,j,k)  = SQRT(N2(i,j,k))
        ENDDO
       ENDDO
      ENDDO

C     N2(k=1) is always zero
      DO j=1-Oly,sNy+Oly
       DO i=1-Olx,sNx+Olx
        N(i,j,1)  = 0.0 _d 0
       ENDDO
      ENDDO

C     Integrate the buoyancy frequency vertically
      DO k=1,Nr
       kp1 = min(k+1,Nr)
       DO j=1-Oly,sNy+Oly
        DO i=1-Olx,sNx+Olx
           BVint(i,j) = BVint(i,j) + hFacC(i,j,k,bi,bj)*drF(k)
     &                         *op5*(MaskC(i,j,k,bi,bj)*N(i,j,k)
     &                         +MaskC(i,j,kp1,bi,bj)*N(i,j,kp1))
        ENDDO
       ENDDO
      ENDDO   

C     Calculate depth-mean buoyancy frequency
      DO j=1-Oly,sNy+Oly
       DO i=1-Olx,sNx+Olx
        IF (r_Low(i,j,bi,bj).NE. 0.0 _d 0) THEN
C         The minus sign is because r_Low<0
          BVzavg(i,j) = -BVint(i,j)/r_Low(i,j,bi,bj)
        ELSE
          BVzavg(i,j) = 0.0 _d 0
        ENDIF
       ENDDO
       ENDDO  

C    Compute slope Burger number and slope-aware modification
      DO j=1-Oly+1,sNy+Oly-1
       DO i=1-Olx,sNx+Olx
        IF (fCori(i,j,bi,bj).NE. 0.0 _d 0) THEN
        SBurger(i,j) = BVzavg(i,j)/fCori(i,j,bi,bj)*abs(dHdy(i,j))
        ELSE
        SBurger(i,j) =  0.0 _d 0  
        ENDIF
        SlopeAwareModi(i,j) = 1/(7.8*SBurger(i,j)**1.5+1)
       ENDDO
      ENDDO

      ENDIF
      
      IF(GMconst_SlopeAware) THEN                    
      DO k=Nr,1,-1
        DO j=1-OLy+1,sNy+OLy-1
          DO i=1-OLx,sNx+OLx
            GEOMK(i,j,k,bi,bj) = GEOM_Const_K*SlopeAwareModi(i,j)
     &                         * maskC(i,j,k,bi,bj)
          ENDDO
        ENDDO
      ENDDO
      ENDIF
      
#ifdef ALLOW_DIAGNOSTICS
      IF (useDiagnostics) THEN
       CALL DIAGNOSTICS_FILL(GEOMK(:,:,1,bi,bj),
     &                       'GEOMkap0', 0,1,0,1,1,myThid)         
      CALL DIAGNOSTICS_FILL(SBurger, 
     &                       'SBurger ', 0,1,0,1,1,myThid)
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
     
      
      
      RETURN
      END