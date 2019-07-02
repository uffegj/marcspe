*  INCLUDE FILE FOR MARCS   version 2005 DJ

* Number of depth points in the model.
      INTEGER   NDP
      PARAMETER ( NDP=100)

* Maximal number of rays for the solution of the transfer equation.
* This should be at least NDP/2 +3.
      INTEGER   NRAYS
      PARAMETER ( NRAYS=80) 

* Equation of state parameters   (for tabs and absko)
* mkomp  = max # of absorption coefficient contributors
* mkompr = max # of temp. indep. absorption coefficient contributors
      integer mkomp,mkompr,ifadim,kfadim
      parameter (mkomp=27,mkompr=16,
     &  ifadim=(mkomp-mkompr)*ndp,kfadim=ifadim*3+mkompr*ndp)

* Number of wavelength points in the  OS calculation
      INTEGER NWL
      PARAMETER (NWL=20000)

* For call to INTERP from OSTABLOOK
      integer mtemp,mpe
      parameter (mtemp=12,mpe=1)

* NUMBER OF ELEMENTS TO DELIVER IONS in JON (and joinabs.dat)
      integer n16,n17
      parameter (n16=16,n17=17)

* Include file for testing Tsuji's molecular equilibrium
* maxmol is maximum number of molecules to be stored in 'takemolec'
* maxmet is maximum number of atoms
* npoint is the number of copmuting temperatures
      integer maxmol,maxmet
      parameter (maxmol=350, maxmet=100)

* maxosmol is the maximum number of molecules for which opacity is included
* (including one number for the atomic lines)
      integer maxosmol
      parameter (maxosmol=30)

* natms = number of atoms considered in Gibbs energy minimalization routines
* nspec = total number of species (atoms, ions, molecules) considered in Gibbs 
*         energy minimalization routines (= 0.5 * number of lines in 
*         /home/ast/uffegj/falkesg/gfits.data
      integer natms,nspec
      parameter(natms=99,nions=127,nspec=896)       ! # elements and species

* Dust
* max_lay = maximum number of DRIFT layers
* max_inc = maximum number of inclusions
* max_eps = maximum number of affected elements
      INTEGER MAX_LAY, MAX_INC, MAX_EPS
      PARAMETER (MAX_LAY = 1000, MAX_INC = 20, MAX_EPS = 10)




