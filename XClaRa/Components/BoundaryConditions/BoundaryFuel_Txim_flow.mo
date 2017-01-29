within Exergy.XClaRa.Components.BoundaryConditions;
model BoundaryFuel_Txim_flow
  "A source defining mass flow, temperature and composition"
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
    powerIn=if massFlowIsLoss then 0 else min(0, fuel_a.m_flow*h_coal),
    powerOut=if massFlowIsLoss then 0 else max(0, fuel_a.m_flow*h_coal),
    powerAux=0) if                                                                                                     contributeToCycleSummary;
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean massFlowIsLoss = true
    "True if mass flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));

  parameter ClaRa.Basics.Media.Fuel.PartialFuel fuelType=simCenter.fuelModel1
    "Coal elemental composition used for combustion"                                          annotation(choices(choice=simCenter.coalModel
        "Coal model 1 as defined in simCenter"),                                              Dialog(group="Fundamental Definitions"));

  parameter Boolean variable_m_flow=false
    "True, if mass flow defined by variable input"                                       annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_T=false
    "True, if temperature defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_xi=false
    "True, if composition defined by variable input"                                      annotation(Dialog(group="Define Variable Boundaries"));

  parameter SI.MassFlowRate m_flow_const=0 "Constant mass flow rate" annotation(Dialog(group="Constant Boundaries", enable= not variable_m_flow));
  parameter SI.Temperature T_const=simCenter.T_amb_start
    "Constant specific temperature of source" annotation(Dialog(group="Constant Boundaries", enable= not hInputIsActive));
  parameter SI.MassFraction xi_const[fuelType.nc-1]=fuelType.defaultComposition
    "Constant composition" annotation(Dialog(group="Constant Boundaries", enable= not variable_xi));
                                                                                            /*zeros(fuelType.nc-1) */
    parameter String LHV_calculationType="predefined" "Calculation type" annotation (
      Dialog(group="Combustion settings"), choices(
      choice="predefined" "Use fixed value for the LHV",
      choice="Verbandsformel" "Calculate the LHV from the Verbandsformel"));

  parameter ClaRa.Basics.Units.EnthalpyMassSpecific LHV_predefined=30e6
    "LHV value for the coal" annotation (Dialog(enable=(LHV_calculationType ==
          "predefined"), group="Combustion settings"));
  parameter Modelica.SIunits.SpecificHeatCapacity cp=fuelType.cp
    "Specific heat capacity of fuel"                                                              annotation (Dialog(group="Combustion settings"));

  outer ClaRa.SimCenter simCenter;
protected
  Modelica.SIunits.MassFlowRate m_flow_in;
  Modelica.SIunits.Temperature T_in;
  Modelica.SIunits.MassFraction xi_in[fuelType.nc-1];
  SI.EnthalpyMassSpecific h_coal;
  //SI.EnthalpyMassSpecific LHV;

public
  ClaRa.Basics.Interfaces.Fuel_outlet fuel_a(final fuelType=fuelType)
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  Modelica.Blocks.Interfaces.RealInput m_flow(value=m_flow_in) if (variable_m_flow)
    "Variable mass flow rate"
    annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
  Modelica.Blocks.Interfaces.RealInput T(value=T_in) if (variable_T)
    "Variable specific temperature"
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
  Modelica.Blocks.Interfaces.RealInput xi[fuelType.nc-1](value=xi_in) if
       (variable_xi) "Variable composition"
    annotation (Placement(transformation(extent={{-120,-80},{-80,-40}})));
equation

  if LHV_calculationType == "predefined" then
    fuel_a.LHV_outflow = LHV_predefined;
    fuel_a.LHV_calculationType=LHV_calculationType;
  elseif LHV_calculationType == "Verbandsformel" then
    fuel_a.LHV_outflow =(33907*fuel_a.xi_outflow[1] + 142324*(fuel_a.xi_outflow[2] - fuel_a.xi_outflow[3]/8.) + 10465*fuel_a.xi_outflow[5] - 2512*((1 - sum(fuel_a.xi_outflow)) + 9*fuel_a.xi_outflow[2]))*1000;
    fuel_a.LHV_calculationType=LHV_calculationType;
  else
    fuel_a.LHV_outflow = LHV_predefined;
    fuel_a.LHV_calculationType=LHV_calculationType;
    assert(fuel_a.LHV_calculationType == "predefined" or fuel_a.LHV_calculationType == "Verbandsformel", "Please check your LHV calculation settings inside boundaries.");
  end if;

  h_coal = actualStream(fuel_a.LHV_outflow) + fuel_a.fuelType.cp*(actualStream(fuel_a.T_outflow) - 298.15);
  fuel_a.cp_outflow=cp;

  if (not variable_m_flow) then
    m_flow_in=m_flow_const;
  end if;
  if (not variable_T) then
    T_in=T_const;
  end if;
  if (not variable_xi) then
    xi_in=xi_const;
  end if;

  fuel_a.T_outflow=T_in;
  fuel_a.m_flow=-m_flow_in;
  fuel_a.xi_outflow=xi_in;

 annotation (Icon(graphics={
        Text(
          extent={{-100,30},{60,-30}},
          lineColor={27,36,42},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="T, xi")}),
                             Diagram(graphics));
end BoundaryFuel_Txim_flow;
