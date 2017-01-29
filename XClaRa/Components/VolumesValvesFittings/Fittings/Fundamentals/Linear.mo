within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals;
model Linear
  extends Fundamentals.BaseDp;
  SI.Pressure dp;
  SI.MassFlowRate m_flow;
  parameter SI.MassFlowRate m_flow_nom=10;
   parameter SI.Pressure dp_nom=1000 "Nominal pressure loss";
equation
  m_flow=m_flow_nom*dp/dp_nom;
end Linear;
