within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
model Quadratic_EN60534_incompressible
  "Quadratic|Kv definition | supercritical flow | incompressible |EN60534"
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.GenericPressureLoss;
  import SI = ClaRa.Basics.Units;
  import SM = ClaRa.Basics.Functions.Stepsmoother;
  parameter Real Kvs(unit="m3/h") = 1
    "|Valve Characteristics|Flow Coefficient at nominal opening (Delta_p_nom = 1e5 Pa, rho_nom=1000 kg/m^3(cold water))";
  Real Kv(unit="m3/h")
    "|Valve Characteristics|Flow Coefficient (Delta_p_nom = 1e5 Pa, rho_nom=1000 kg/m^3(cold water))";
  parameter Real F_L= 0.75
    "|Valve Characteristics|Relative Pressure drop without fittings (see docu for typical values)";
  parameter SI.MassFlowRate m_flow_nominal= Kvs/1000/3600
    "|Valve Characteristics|Only for homotopy-based initialisation: Nominal mass flowrate at full opening";
  parameter Real diameter_inlet(unit="mm") = 10
    "|Valve Characteristics|Inlet fitting's diameter";
  parameter Real diameter_valve(unit="mm") = 10
    "|Valve Characteristics|Valve diameter";
  parameter Real diameter_outlet(unit="mm") = 10
    "|Valve Characteristics|Outlet fitting's diameter";

  parameter SI.Pressure Delta_p_eps= 100
    "|Expert Settings||Small pressure difference for linearisation around zeor flow";
  SI.Pressure Delta_p_choke;
//  Real Delta_p_ "(p_inlet-p_outlet)/p_inlet";

protected
  constant Real K1= 0.0016 "Coeffitient for calculation of F_P";

  Real F_P = 1/sqrt(1+(((1-diameter_valve^4/diameter_inlet^4) - (1-diameter_valve^4/diameter_outlet^4) + 0.5*(1-(diameter_valve/diameter_inlet)^2)^2 + (1-(diameter_valve/diameter_outlet)^2)^2)*Kv^2/diameter_valve^4)/K1)
    "Pipe geometry correction factor";
  SI.Pressure p_in = max(iCom.p_in, iCom.p_out);
  Real F_LP = F_L / sqrt(1+(F_L^2*((1-diameter_valve^4/diameter_inlet^4) + 0.5*(1-(diameter_valve/diameter_inlet)^2)^2)/K1*(Kv^2/diameter_valve^4)));
equation
  gamma = 2e30 "is not used";
 //____________ Pressure drop in design flow direction_______________
//  Delta_p_ = if checkValve then noEvent(max(0, abs(Delta_p))/p_in) else abs(Delta_p)/p_in;
  Delta_p_choke =(F_LP/F_P)^2*(p_in - (0.96-0.28*sqrt(0.0737/221.2))*0.0737);
//____________Hydraulics____________________________________________
  Kv = aperture_ * Kvs;

  //m_flow =  noEvent(if checkValve then noEvent(max(0,noEvent(if Delta_p >1 then (if useHomotopy then homotopy(1/3600*0.1*Kv*F_P*ClaRa.Basics.Functions.ThermoRoot(min(Delta_p,Delta_p_choke), Delta_p_eps)*sqrt(iCom.rho_in), aperture_*m_flow_nominal)
  //                                     else 1/3600*0.1*Kv*F_P*
  //  ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in)) else 0))) else if useHomotopy then homotopy(1/3600*0.1*Kv*F_P*ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in), aperture_*m_flow_nominal)
  //                                     else 1/3600*0.1*Kv*F_P*
  //  ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in));
      m_flow =  if checkValve then if useHomotopy then homotopy(SM(Delta_p_eps/10,0, Delta_p) * 1/3600*0.1*Kv*F_P*ClaRa.Basics.Functions.ThermoRoot(min(Delta_p,Delta_p_choke), Delta_p_eps)*sqrt(iCom.rho_in), aperture_*m_flow_nominal)
                                       else SM(Delta_p_eps/10,0, Delta_p)*1/3600*0.1*Kv*F_P*
    ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in) else if useHomotopy then homotopy(1/3600*0.1*Kv*F_P*ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in), aperture_*m_flow_nominal)
                                       else 1/3600*0.1*Kv*F_P*
    ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in);

end Quadratic_EN60534_incompressible;
