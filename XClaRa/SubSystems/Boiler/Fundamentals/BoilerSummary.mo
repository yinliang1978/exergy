within Exergy.XClaRa.SubSystems.Boiler.Fundamentals;
record BoilerSummary
  "A summary for boilers - takes only the water steam side at the connectors into account"
  extends ClaRa.Basics.Icons.RecordIcon;

  input Modelica.SIunits.Pressure p_feed "Feedwater inlet pressure" annotation(Dialog(group="Feedwater"));
  input Modelica.SIunits.SpecificEnthalpy h_feed
    "Feedwater inlet specific enthalpy"
                                       annotation(Dialog(group="Feedwater"));
  input Modelica.SIunits.MassFlowRate m_flow_feed "Feedwater inlet mass flow"
                                                                             annotation(Dialog(group="Feedwater"));

  input Modelica.SIunits.Pressure p_LS "Live steam pressure"
                                                            annotation(Dialog(group="Live steam"));
  input Modelica.SIunits.SpecificEnthalpy h_LS "Live steam specific enthalpy"
                                                                             annotation(Dialog(group="Live steam"));
  input Modelica.SIunits.MassFlowRate m_flow_LS "Live steam mass flow"
                                                                      annotation(Dialog(group="Live steam"));

  input Modelica.SIunits.Pressure p_cRH "Cold reheat pressure"
                                                              annotation(Dialog(group="Cold Reheat"));
  input Modelica.SIunits.SpecificEnthalpy h_cRH "Cold reheat specific enthalpy"
                                                                                annotation(Dialog(group="Cold Reheat"));
  input Modelica.SIunits.MassFlowRate m_flow_cRH "Cold reheat mass flow"
                                                                        annotation(Dialog(group="Cold Reheat"));

  input Modelica.SIunits.Pressure p_hRH "Hot reheat pressure"
                                                             annotation(Dialog(group="Hot Reheat"));
  input Modelica.SIunits.SpecificEnthalpy h_hRH "Hot reheat specific enthalpy"
                                                                              annotation(Dialog(group="Hot Reheat"));
  input Modelica.SIunits.MassFlowRate m_flow_hRH "Hot reheat mass flow"
                                                                       annotation(Dialog(group="Hot Reheat"));

end BoilerSummary;
