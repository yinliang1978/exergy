within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
model LinearNominalPoint "Linear|Nominal operation point | subcritical flow"
//   "A linear pressure loss using a constant pressure loss coefficient"
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.GenericPressureLoss;
  import SI = ClaRa.Basics.Units;
  import SM = ClaRa.Basics.Functions.Stepsmoother;
  parameter SI.PressureDifference Delta_p_nom = 1e5
    "|Valve Characteristics|Nominal pressure difference for Kv definition";

  parameter SI.MassFlowRate m_flow_nom= 1
    "|Valve Characteristics|Nominal mass flow rate";
 // parameter SI.PressureDifference Delta_p_check=0;

 // parameter SI.PressureDifference Delta_p_hyst=0;

 // Boolean ValveOpen;

equation
  gamma = 2e30 "is not used";
//  m_flow = noEvent(if checkValve then if Delta_p >0 then aperture_* m_flow_nom * Delta_p/Delta_p_nom else 0 else aperture_* m_flow_nom * Delta_p/Delta_p_nom);

//   m_flow = noEvent(if checkValve then if Delta_p >=0 then aperture_* m_flow_nom * Delta_p/Delta_p_nom else (if iCom.opening_leak_ >=1e-6 then iCom.opening_leak_* m_flow_nom * Delta_p/Delta_p_nom else Modelica.Constants.eps) else
//                                                                                                     aperture_* m_flow_nom * Delta_p/Delta_p_nom);

//ValveOpen= if m_flow>0 then true else false;

// m_flow = noEvent(
// if checkValve then
//   if Delta_p >Delta_p_check + Delta_p_hyst and not pre(ValveOpen) or Delta_p >= Delta_p_check and pre(ValveOpen)  then
//   aperture_* m_flow_nom * (Delta_p-Delta_p_check)/Delta_p_nom
//   else
//   iCom.opening_leak_* m_flow_nom * Delta_p/Delta_p_nom
//  else
// aperture_* m_flow_nom * (Delta_p-Delta_p_check)/Delta_p_nom);

//m_flow= if checkValve then SM(50,0, Delta_p)*aperture_* m_flow_nom * (Delta_p)/Delta_p_nom + SM(-10,0, Delta_p)*(-1e-5) else aperture_* m_flow_nom * (Delta_p)/Delta_p_nom;
m_flow= if checkValve then SM(50,0, Delta_p)*aperture_* m_flow_nom * (Delta_p)/Delta_p_nom else aperture_* m_flow_nom * (Delta_p)/Delta_p_nom;
//or
end LinearNominalPoint;
