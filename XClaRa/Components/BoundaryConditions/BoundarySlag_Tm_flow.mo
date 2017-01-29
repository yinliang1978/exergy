within Exergy.XClaRa.Components.BoundaryConditions;
model BoundarySlag_Tm_flow "A source defining mass flow and temperature"
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                        //
//                                                                           //
// Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
// Copyright © 2013-2016, DYNCAP/DYNSTART research team.                     //
//___________________________________________________________________________//
// DYNCAP and DYNSTART are research projects supported by the German Federal //
// Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
// The research team consists of the following project partners:             //
// Institute of Energy Systems (Hamburg University of Technology),           //
// Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
// TLK-Thermo GmbH (Braunschweig, Germany),                                  //
// XRG Simulation GmbH (Hamburg, Germany).                                   //
//___________________________________________________________________________//

  extends ClaRa.Basics.Icons.FlowSource;
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=if massFlowIsLoss then 0 else min(0, slag_outlet.m_flow*h_slag),
    powerOut=if massFlowIsLoss then 0 else max(0, slag_outlet.m_flow*h_slag),
    powerAux=0) if                                                                                                     contributeToCycleSummary;
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean massFlowIsLoss = true
    "True if mass flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));

  parameter ClaRa.Basics.Media.Fuel.PartialSlag slagType=simCenter.slagModel
    "Medium to be used"                                                                                               annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  parameter Boolean variable_m_flow=false
    "True, if mass flow defined by variable input"                                       annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_T=false
    "True, if temperature defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
  parameter SI.MassFlowRate m_flow_const=0 "Constant mass flow rate" annotation(Dialog(group="Constant Boundaries", enable= not variable_m_flow));
  parameter SI.Temperature T_const=simCenter.T_amb_start
    "Constant specific temperature of source" annotation(Dialog(group="Constant Boundaries", enable= not hInputIsActive));

  outer ClaRa.SimCenter simCenter;
protected
  Modelica.SIunits.MassFlowRate m_flow_in;
  Modelica.SIunits.Temperature T_in;
  SI.EnthalpyMassSpecific h_slag;

public
        ClaRa.Basics.Interfaces.Slag_outlet slag_outlet(final slagType=slagType)
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  Modelica.Blocks.Interfaces.RealInput m_flow(value=m_flow_in) if (variable_m_flow)
    "Variable mass flow rate"
    annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
  Modelica.Blocks.Interfaces.RealInput T(value=T_in) if (variable_T)
    "Variable specific temperature"
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));

equation
  h_slag = slag_outlet.m_flow*slagType.cp*(actualStream(slag_outlet.T_outflow) - 298.15);

  if (not variable_m_flow) then
    m_flow_in=m_flow_const;
  end if;
  if (not variable_T) then
    T_in=T_const;
  end if;

  slag_outlet.T_outflow=T_in;
  slag_outlet.m_flow=-m_flow_in;

 annotation (Icon(graphics={
        Text(
          extent={{-100,30},{60,-30}},
          lineColor={27,36,42},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="T, xi")}),
                             Diagram(graphics));
end BoundarySlag_Tm_flow;
