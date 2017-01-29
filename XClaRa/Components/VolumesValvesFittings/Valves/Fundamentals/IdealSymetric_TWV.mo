within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
model IdealSymetric_TWV
  "| No Pressure Dependency | Opening Characteristics | Symetrical |"
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.Basic_TWV;
  extends
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.TWV_L2;

  import SI = ClaRa.Basics.Units;

  SI.MassFlowRate m_flow_1 "Mass flow through straight flange";
  SI.MassFlowRate m_flow_2 "Mass flow through shunt flange";
  SI.Pressure p_out2;

//  SI.PressureDifference Delta_p "Pressure difference p_in - p_out";
//  Real gamma;

equation
  m_flow_1 = iCom.m_flow_in*aperture_;
  m_flow_2 = iCom.m_flow_in*(1-aperture_);
  p_out2 = iCom.p_in;

end IdealSymetric_TWV;
