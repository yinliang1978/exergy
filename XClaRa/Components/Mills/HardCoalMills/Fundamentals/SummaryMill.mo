within Exergy.XClaRa.Components.Mills.HardCoalMills.Fundamentals;
model SummaryMill "A summary record for a roller bowl mill"

extends ClaRa.Basics.Icons.RecordIcon;
  input ClaRa.Basics.Units.MassFlowRate m_flow_coal_in
    "Coal mass flow entering the mill"                                                    annotation(Dialog);
  input ClaRa.Basics.Units.MassFlowRate m_flow_coal_out
    "Coal mass flow leaving the mill"                                                     annotation(Dialog);
  input ClaRa.Basics.Units.MassFlowRate m_flow_air_in
    "Coal mass flow entering the mill"                                                    annotation(Dialog);
  input ClaRa.Basics.Units.MassFlowRate m_flow_air_out
    "Coal mass flow leaving the mill"                                                    annotation(Dialog);
  input ClaRa.Basics.Units.MassFlowRate m_flow_tot_in
    "Total mass flow entering the mill"                                                    annotation(Dialog);
  input ClaRa.Basics.Units.MassFlowRate m_flow_tot_out
    "Total mass flow leaving the mill"                                                    annotation(Dialog);

  input ClaRa.Basics.Units.Temperature T_out "Classifier temperature" annotation(Dialog);
  input ClaRa.Basics.Units.Temperature T_coal_in "Coal inlet temperature" annotation(Dialog);
  input ClaRa.Basics.Units.Temperature T_air_in "Primary air inlet temperature"
                                                                                annotation(Dialog);
  input ClaRa.Basics.Units.RPM rpm_classifier "Classifier speed" annotation(Dialog);
  input ClaRa.Basics.Units.Power P_grind "Power consumed for grinding" annotation(Dialog);
  input ClaRa.Basics.Units.Mass mass_coal "Total coal mass in the mill" annotation(Dialog);

  input ClaRa.Basics.Units.MassFraction xi_coal_h2o_in
    "Water mass fraction of incoming coal"                                                    annotation(Dialog);
  input ClaRa.Basics.Units.MassFraction xi_coal_h2o_out
    "Water mass fraction of leaving coal"                                                     annotation(Dialog);

  input ClaRa.Basics.Units.MassFraction xi_air_h2o_in
    "Water mass fraction of incoming air"                                                   annotation(Dialog);
  input ClaRa.Basics.Units.MassFraction xi_air_h2o_out
    "Water mass fraction of leaving air"                                                    annotation(Dialog);
  input ClaRa.Basics.Units.MassFraction xi_air_h2o_sat
    "Max. water mass fraction at air path"                                                    annotation(Dialog);

end SummaryMill;
