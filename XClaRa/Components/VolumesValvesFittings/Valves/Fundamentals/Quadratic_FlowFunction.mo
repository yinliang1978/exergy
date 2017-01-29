within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
model Quadratic_FlowFunction "Quadratic|zeta definition | supercritical flow"
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.GenericPressureLoss;
  import SI = ClaRa.Basics.Units;
  import SM = ClaRa.Basics.Functions.Stepsmoother;
  outer Boolean checkValve;
  outer Boolean useHomotopy;

  parameter SI.Area A_cross= 1
    "|Valve Characteristics|Valve inlet cross section";
  parameter Real zeta= 100 "|Valve Characteristics|Pressure Loss coefficient";
  final parameter Real Kvs(unit="m3/h") = A_cross*sqrt(2/zeta) *sqrt(1e5/1000)*3600
    "|Valve Characteristics|Flow Coefficient at nominal opening (Delta_p_nom = 1e5 Pa, rho_nom=1000 kg/m^3(cold water))";

  parameter SI.MassFlowRate m_flow_nominal= Kvs/1000/3600
    "|Valve Characteristics|Only for homotopy-based initialisation: Nominal mass flowrate at full aperture_";

  parameter SI.Pressure Delta_p_eps= 100
    "|Expert Settings||Small pressure difference for linearisation around zeor flow";

//  Real Delta_p_ "(p_inlet-p_outlet)/p_inlet";

protected
  Real Psi;
 SI.Pressure p_in = max(iCom.p_in, iCom.p_out);
equation
 gamma = if (checkValve == true and iCom.opening_leak_ <= 0) or iCom.opening_ < iCom.opening_leak_ then iCom.gamma_in else (if useHomotopy then homotopy(ClaRa.Basics.Functions.Stepsmoother(10, -10, Delta_p)*iCom.gamma_in + ClaRa.Basics.Functions.Stepsmoother(-10, 10, Delta_p)*iCom.gamma_out, iCom.gamma_in) else ClaRa.Basics.Functions.Stepsmoother(10, -10, Delta_p)*iCom.gamma_in + ClaRa.Basics.Functions.Stepsmoother(-10, 10, Delta_p)*iCom.gamma_out);

//____________Hydraulics____________________________________________

  Psi = noEvent(if min(iCom.p_in/iCom.p_out,iCom.p_out/iCom.p_in)> (2/(gamma+1))^(gamma/(gamma-1)) then sqrt(gamma/(gamma-1))*sqrt(min(iCom.p_in/iCom.p_out, iCom.p_out/iCom.p_in)^(2/gamma) - min(iCom.p_in/iCom.p_out,iCom.p_out/iCom.p_in)^((gamma+1)/gamma)) else sqrt(gamma/(gamma+1))*(2/(gamma+1))^(1/(gamma-1)));
 // m_flow = noEvent(if checkValve then if Delta_p >0 then if useHomotopy then homotopy(sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps), aperture_*Psi*m_flow_nominal)
 //else
 //    sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps) else 0 else if useHomotopy then homotopy(sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps), aperture_*Psi*m_flow_nominal)
 //else
 //    sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps));
       m_flow = if checkValve then if useHomotopy then homotopy(SM(Delta_p_eps/10,0, Delta_p)*sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps), aperture_*Psi*m_flow_nominal)
 else
     SM(Delta_p_eps/10,0, Delta_p)*sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps) else if useHomotopy then homotopy(sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps), aperture_*Psi*m_flow_nominal)
 else
     sign(iCom.p_in-iCom.p_out)*aperture_*Psi *ClaRa.Basics.Functions.ThermoRoot(2*iCom.rho_in*p_in/zeta, 2*iCom.rho_in/zeta* Delta_p_eps);

end Quadratic_FlowFunction;
