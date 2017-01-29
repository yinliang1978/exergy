within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals;
model NoFriction
  extends Fundamentals.BaseDp;
  SI.Pressure dp;
  SI.MassFlowRate m_flow;
equation
  dp=0;
end NoFriction;
