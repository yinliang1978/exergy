within Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals;
partial model Basic_TWV "Basic pressure loss model for three way valves"
  extends ClaRa.Basics.Icons.Delta_p;
  import SI = ClaRa.Basics.Units;

  parameter Real CharLineValve_[:, :]=[0, 0; 1, 1]
    "|Valve Characteristics|Effective apperture as function of valve position in p.u.";
  parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
    "|Expert Settings|Numerical Robustness|Smoothness of table interpolation";

  outer
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.ICom_TWV
    iCom;
//   outer Boolean checkValve;
//   outer Boolean useHomotopy;

  SI.MassFlowRate m_flow_1 "Mass flow through straight flange";
  SI.MassFlowRate m_flow_2 "Mass flow through shunt flange";
  Real aperture_ "Effective apperture";

protected
  Exergy.XClaRa.Components.Utilities.Blocks.ParameterizableTable1D table(
    table=CharLineValve_,
    columns={2},
    smoothness=smoothness) "External table object";

equation
  table.u[1] = iCom.opening_;
  aperture_ = min(1, max(0, table.y[1]));

end Basic_TWV;
