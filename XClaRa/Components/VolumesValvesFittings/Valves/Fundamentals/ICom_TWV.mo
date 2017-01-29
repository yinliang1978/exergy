within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
record ICom_TWV
  extends ClaRa.Basics.Icons.RecordIcon;
  import SI = ClaRa.Basics.Units;

  SI.MassFlowRate m_flow_in "Inlet mass flow rate" annotation(Dialog);
  SI.Pressure p_in "Inlet pressure" annotation(Dialog);
  SI.Pressure p_out1 "Outlet pressure at straight flange" annotation(Dialog);
  SI.Pressure p_out2 "Outlet pressure at shunt flange" annotation(Dialog);
  SI.DensityMassSpecific rho_in "Outlet density at straight flange" annotation(Dialog);
  SI.DensityMassSpecific rho_out1 "Outlet density at straight flange" annotation(Dialog);
  SI.DensityMassSpecific rho_out2 "Outlet density at shunt flange" annotation(Dialog);
  Real opening_ "Opening of the valve (p.u.)" annotation(Dialog);

end ICom_TWV;
