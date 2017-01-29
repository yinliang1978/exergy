within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
model QuadraticFrictionNominalPointSymetric_TWV
  "| Quadratic Pressure Dependency | Nominal Point | Opening Characteristics | Symetrical |"
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.Basic_TWV;
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.TWV_L1;
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.TWV_L2;

  import SI = ClaRa.Basics.Units;

  parameter SI.Area effectiveFlowArea1=7.85e-3
    "Effective flow area for straight outlet"
    annotation(Dialog(group="Valve Characteristics"));

  parameter SI.Pressure Delta_p_nom[2] = {1e5,1e5}
    "|Valve Characteristics|Nominal pressure difference for Kv definition";
  parameter SI.MassFlowRate m_flow_nom[2] = {1,1}
    "|Valve Characteristics|Nominal mass flow rate";
  parameter Boolean useStabilisedMassFlow=false
    "|Expert Settings|Numerical Robustness|";
  parameter SI.Time Tau= 0.001 "Time Constant of Stabilisation" annotation(Dialog(tab="Expert Settings", group = "Numerical Robustness", enable=useStabilisedMassFlow));
  parameter SI.PressureDifference Delta_p_smooth = 100
    "Below this value, root function is approximated linearly"                                                    annotation(Dialog(tab = "Expert Settings", group="Numerical Robustness"));

  final parameter Real Kvs[2](unit="m3/h") = 3600 * m_flow_nom./rho_in_nom .* sqrt(rho_in_nom/1000*1e5./Delta_p_nom)
    "|Valve Characteristics|Flow Coefficient at nominal opening (Delta_p_nom = 1e5 Pa, rho_nom=1000 kg/m^3(cold water))";
  final parameter SI.DensityMassSpecific rho_in_nom= 1000
    "|Valve Characteristics|Nominal density for Kv definition";

  Real Kv[2](unit="m3/h")
    "|Valve Characteristics|Flow Coefficient (Delta_p_nom = 1e5 Pa, rho_nom=1000 kg/m^3(cold water))";
  SI.Pressure Delta_p[2](start={10,10}) "Pressure differences";
equation
//////////// Simple hydraulics: ///////////////////////////////
  if useStabilisedMassFlow==false then
    Delta_p[1] = iCom.p_in - iCom.p_out1;
    Delta_p[2] = iCom.p_in - iCom.p_out2;
  else
    der(Delta_p[1]) = (iCom.p_in - iCom.p_out1 - Delta_p[1])/Tau;
    der(Delta_p[2]) = (iCom.p_in - iCom.p_out2 - Delta_p[1])/Tau;
  end if;
   Kv[1] = aperture_ * Kvs[1];
   Kv[2] = (1-aperture_) * Kvs[2];
   m_flow_1 = Kv[1]/3600 * rho_in_nom * ClaRa.Basics.Functions.ThermoRoot(Delta_p[1]/1e5, Delta_p_smooth/1e5)*sqrt(1000/rho_in_nom);
   m_flow_2 = Kv[2]/3600 * rho_in_nom * ClaRa.Basics.Functions.ThermoRoot(Delta_p[2]/1e5, Delta_p_smooth/1e5)*sqrt(1000/rho_in_nom);

initial equation
  if useStabilisedMassFlow then
    Delta_p[1] = iCom.p_in - iCom.p_out1;
    Delta_p[2] = iCom.p_in - iCom.p_out2;
  end if;
end QuadraticFrictionNominalPointSymetric_TWV;
