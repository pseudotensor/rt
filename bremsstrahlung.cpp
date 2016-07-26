// THERMAL BREMSSTRAHLUNG
// Roman Gold, started Jun 28th 2016

// REF: http://arxiv.org/pdf/1606.08192.pdf section 4.2 eqs (26)-(34)

doub F_ei(doub Theta_e) 
  {
    if (Theta_e<1) return 4.* sqrt(2.*Theta_e/PI/PI/PI) * ( 1. + 1.781*pow(Theta_e,1.34) );
    else         return 9.*Theta_e/2./PI*(1.5 * log(1.123*Theta_e + 0.48) );
  }

doub gbar(doub x)
{ // auxiliary fct for Bremsstrahlug emissivity
  if (x<1) return sqrt(3)/PI*log(2.246/x);
  else     return sqrt(3./PI/x);
}


doub j_br(doub rho_cgs, doub t_e, doub nu){

  // COMPUTE BREMSSTRAHLUNG EMISSIVITY (IN CGS UNITS)
  // e:electron i:ion assume fully ionized plasma of electrons and protons

  printf(YELLOW"[bremsstrahlung.cpp]: "RESET"rho_cgs=%f t_e=%f nu=%f\n",rho_cgs, t_e, nu);

  doub x = h_planck * nu/kb/t_e; // eq (32)
  doub Theta_e = kb*t_e/me/cc/cc; // eq (28) dim less electron temperature 

  doub C1=7.028e25;
  doub C2=9.334e25;

  doub q_ee = 0.; // eq (30) units: [erg/cm^3/s]
  if (Theta_e < 1) q_ee = C1 * rho_cgs*rho_cgs * pow(Theta_e,1.5) * (1.+1.1*Theta_e + Theta_e*Theta_e - 1.25*pow(Theta_e,2.5)); 
  else q_ee = C2 * rho_cgs*rho_cgs * Theta_e * (log(1.123*Theta_e + 1.28));

  doub q_ei = 3.013e25 * rho_cgs*rho_cgs * F_ei(Theta_e); // eq (27) [erg/cm^3/s]
  if (Theta_e < 1) q_ei = C1 * rho_cgs*rho_cgs * pow(Theta_e,1.5) * (1.+1.1*Theta_e + Theta_e*Theta_e - 1.25*pow(Theta_e,2.5)); 
  else q_ei = C2 * rho_cgs*rho_cgs * Theta_e * (log(1.123*Theta_e + 1.28));

  doub q_br = q_ei + q_ee; // eq (26)
  doub j_br = 1./4./PI/nu * q_br * x * exp(-x) * gbar(x); // eq (31) bremsstrahlung emissivity

  return j_br;
}
