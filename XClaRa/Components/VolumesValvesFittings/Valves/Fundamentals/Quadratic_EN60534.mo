within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
model Quadratic_EN60534
  "Quadratic|Kv definition | supercritical flow | compressible |EN60534"
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.GenericPressureLoss;
  import SI = ClaRa.Basics.Units;
  import SM = ClaRa.Basics.Functions.Stepsmoother;
  parameter Real Kvs(unit="m3/h") = 1
    "|Valve Characteristics|Flow Coefficient at nominal opening (Delta_p_nom = 1e5 Pa, rho_nom=1000 kg/m^3(cold water))";
  Real Kv(unit="m3/h")
    "|Valve Characteristics|Flow Coefficient (Delta_p_nom = 1e5 Pa, rho_nom=1000 kg/m^3(cold water))";
  parameter Real x_T= 0.75
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

  Real Delta_p_ "(p_inlet-p_outlet)/p_inlet";

protected
  constant Real K1= 0.0016 "Coefficient for calculation of F_P";
  constant Real K2= 0.0018 "Coefficient for calculation of X_TP";
  Real Y "Expansion factor";
   Real F_P = 1/sqrt(1+(((1-diameter_valve^4/diameter_inlet^4) - (1-diameter_valve^4/diameter_outlet^4) + 0.5*(1-(diameter_valve/diameter_inlet)^2)^2 + (1-(diameter_valve/diameter_outlet)^2)^2)*Kv^2/diameter_valve^4)/K1)
    "Pipe geometry correction factor";
  SI.Pressure p_in = max(iCom.p_in, iCom.p_out);
  Real x_TP = (x_T/F_P^2) / (1+(x_T*((1-diameter_valve^4/diameter_inlet^4) - (1-diameter_valve^4/diameter_outlet^4) + 0.5*(1-(diameter_valve/diameter_inlet)^2)^2 + (1-(diameter_valve/diameter_outlet)^2)^2)/K2*(Kv^2/diameter_valve^4)));
equation
 gamma = if (checkValve == true and iCom.opening_leak_ <= 0) or iCom.opening_ < iCom.opening_leak_ then iCom.gamma_in else (if useHomotopy then homotopy(ClaRa.Basics.Functions.Stepsmoother(10, -10, Delta_p)*iCom.gamma_in + ClaRa.Basics.Functions.Stepsmoother(-10, 10, Delta_p)*iCom.gamma_out, iCom.gamma_in) else ClaRa.Basics.Functions.Stepsmoother(10, -10, Delta_p)*iCom.gamma_in + ClaRa.Basics.Functions.Stepsmoother(-10, 10, Delta_p)*iCom.gamma_out);

 //____________ Pressure drop in design flow direction_______________
  Delta_p_ = if checkValve then noEvent(max(0, Delta_p)/p_in) else abs(Delta_p)/p_in;

//____________Hydraulics____________________________________________
  Kv = aperture_ * Kvs;
  //Y= noEvent(min(1,max(2/3, 1-Delta_p_/(3*iCom.gamma/1.4*(x_T/F_P)))));
  Y= noEvent(min(1,max(2/3, 1-Delta_p_/(3*gamma/1.4*x_TP))));
 // m_flow =  noEvent(if checkValve then if Delta_p >0 then if useHomotopy then homotopy(if aperture_ > 1e-6 then 1/3600*0.1*Kv*F_P*Y*ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in) else 0, aperture_*m_flow_nominal)
 //                                      else 1/3600*0.1*Kv*F_P*Y*
 //   ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in) else 0 else if useHomotopy then homotopy(1/3600*0.1*Kv*F_P*Y*ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in), aperture_*m_flow_nominal)
 //                                      else 1/3600*0.1*Kv*F_P*Y*
 //   ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in));
      m_flow = if checkValve then if useHomotopy then homotopy(SM(Delta_p_eps/10,0, Delta_p) *1/3600*0.1*Kv*F_P*Y*ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in), aperture_*m_flow_nominal)
 else
     SM(Delta_p_eps/10,0, Delta_p) * 1/3600*0.1*Kv*F_P*Y*
    ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in) else if useHomotopy then homotopy(1/3600*0.1*Kv*F_P*Y*ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in), aperture_*m_flow_nominal)
                                       else 1/3600*0.1*Kv*F_P*Y*
    ClaRa.Basics.Functions.ThermoRoot(Delta_p, Delta_p_eps)*sqrt(iCom.rho_in);

end Quadratic_EN60534;
