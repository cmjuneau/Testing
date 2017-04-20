!  ========================================================

  subroutine setheader
    implicit none
    
!
!  This subroutine will create a list of section headers
!  used in CEM.  These headers will indicate to the code
!  when to read output values and store appropriately.

    integer :: headmult
    
    character(LEN = 300) :: header, mischeader
    
    logical :: makeheader, makemisc
    
    
    common  /header/ header(100), mischeader(50), headmult
! Initialize several variables
		header = "- - - - - NO HEADER INFORMATION - - - - -"
    mischeader = "- - - - - NO MISC. INFORMATION - - - - -"
    makeheader = .true.
    makemisc = .true.
    
! Header information
    if(makeheader) then
      header(1)  = " Part.     Multiplicities           Yields [mb]     <TKE> [MeV]"
      header(2)  = " Yields of different channels (with > 1 mb):"
      header(3)  = " Yields of several summed channels:"
      header(4)  = " *************** Nuclide yields [mb]  (zero values suppressed) *****************"
      header(5)  = "Mass yield [mb] and the mean and variance of the kinetic energy [MeV]"
      header(6)  = "Charge yield [mb] and the mean and variance of the  kinetic energy [MeV]"
      header(7)  = " --------------------------- Energy Spectrum [mb/MeV] --------------------------"
      header(8)  = " ----------------Normalized Energy Probability Spectrum [1/MeV] ----------------"
      header(9)  = " ------------------------ Angular Distributions [mb/sr] ------------------------"
      header(10) = " Double differential cross sections [mb/MeV/sr];" ! Lab angle proceeds
      header(11) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    n"
      header(12) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    p"
      header(13) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    d"
      header(14) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of    t"
      header(15) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  He3"
      header(16) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  He4"
      header(17) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  He6"
      header(18) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li6"
      header(19) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li7"
      header(20) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li8"
      header(21) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Li9"
      header(22) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Be7"
      header(23) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  Be9"
      header(24) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of Be10"
      header(25) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  B10"
      header(26) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  B11"
      header(27) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  C11"
      header(28) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of  C12"
      header(29) = "          Angular distribution of produced fragments dS/dOm [mb/sr] for energy range(MeV)"
      header(30) = "                         Double differential cross-section d2S/dTdO (mb/MeV/sr) of"
      header(31) = "          Spectra ( dS/dT/dO ) of produced particles at theta="
      header(32) = "          Angle integrated spectra of produced particles"
      header(33) = "          Energy integ. ang. distrib. of prod. particles"
      header(34) = "          Energy integ. ang. distrib. of prod. particles"
      header(35) = "          Angular distribution of produced fragments dS/dOm [mb/sr] for energy range(MeV)"
      header(36) = " Int. x sec "
      header(37) = " Int. xsec "
    endif
    
! Misc. / Complimentary Header
    if(makemisc) then
      mischeader(1)  = " *****************************************************************"
      mischeader(2)  = " End of nuclide yields."
      mischeader(3)  = " of residual nuclei:" ! used with Mass and Charge yields
! For spectra
      mischeader(4)  = "+/-" ! for uncertainties
      mischeader(5)  = "  [deg.]"
      mischeader(6)  = " Lab. angle =" ! for double differential cross sections
      mischeader(7)  = " **********************************  neutrons  *********************************"
      mischeader(8)  = " **********************************  protons   *********************************"
      mischeader(9)  = " ********************************** deuterons  *********************************"
      mischeader(10) = " **********************************  tritons   *********************************"
      mischeader(11) = " **********************************  Helium-3  *********************************"
      mischeader(12) = " **********************************   alphas   *********************************"
      mischeader(13) = " ********************************** neg. pions *********************************"
      mischeader(14) = " ********************************** neut pions *********************************"
      mischeader(15) = " ********************************** pos. pions *********************************"
      mischeader(16) = "     Energy spectrum from   0.0 to  396.0 MeV (zero values suppressed)."
      mischeader(17) = "    Tn  [MeV]            Total                Cascade               Precompound        Total Evaporation"
      mischeader(18) = "    Tp  [MeV]            Total                Casc"
      mischeader(19) = "    Tt  [MeV]           Total               Coalescence            Precompound         Total Evaporation"
      mischeader(20) = "    THe3[MeV]           Total               Coalescence            Precompound         Total Evaporation"
      mischeader(21) = "    THe4[MeV]           Total               Coalescence            Precompound         Total Evaporation"
      mischeader(22) = "    Tpi-[MeV]      Total = Cascade"
      mischeader(23) = "    Tpi0[MeV]      Total = Cascade"
      mischeader(24) = "    Tpi+[MeV]      Total = Cascade"
      mischeader(25) = "  Ang.n           Total                Cascade               Precompound         Total Evaporation"
      mischeader(26) = "  Ang.p           Total                Cascade               Precompound         Total Evaporation"
      mischeader(27) = "  Ang.t           Total              Coalescence             Precompound         Total Evaporation"
      mischeader(28) = "  Ang.He3         Total              Coalescence             Precompound         Total Evaporation"
      mischeader(29) = "  Ang.He4         Total              Coalescence             Precompound         Total Evaporation"
      mischeader(30) = "  Ang.pi-[deg.]     Total = Cascade"
      mischeader(31) = "  Ang.pi0[deg.]     Total = Cascade"
      mischeader(32) = "  Ang.pi+[deg.]     Total = Cascade"
!      mischeader(33) = "  T(MeV)/angle:  12.        15.        20.        35.        50.        65.        80.        90.       100.       120.    dS/dT(mb/MeV)"
      mischeader(34) = "            prod. xsec for all energies="
!      mischeader(35) = "Tmin-Tmax   0.- 160.   2.- 160.   3.- 215.   3.- 250.   2.- 580.   2.- 650.   2.-2500.   4.-2500.   4.-2500.   4.-2500."
!      mischeader(36) = "ang1-ang2      n          p          d          t        He3        He4        He6        Li6        Li7        Li8"
!      mischeader(37) = "ang1-ang2    Li9        Be7        Be9       Be10         B9        B10        B11        B12        C11        C12"
!      mischeader(38) = "ang1-ang2    C13        C14        Z=7        Z=8        Z=9       Z=10       Z=11       Z=12       Z=13       Z=14"
      mischeader(39) = " The program called Fermi breakup"
      mischeader(40) = "- - - - NO INFORMATION - - - -"
      mischeader(41) = "*******"
      mischeader(42) = " energ. int."
      mischeader(43) = "  T(GeV)      He6        Li6        Li7        Be7        Be9        Be10         B8         B10        B11"
      mischeader(44) = "              He6        Li6        Li7        Be7        Be9        Be10         B8         B10        B11"
      mischeader(45) = " theta      He6        Li6        Li7        Be7        Be9        Be10         B8         B10        B11"
      mischeader(46) = "  T(GeV)    PI-MINUS    PI-PLUS   NEUTRONS    PROTONS   DEUTRONS    TRITONS       HE-3       HE-4"
    endif
    
    return
  end subroutine setheader

!  ========================================================