within Exergy.XClaRa.Components.BoundaryConditions;
model BoundaryGas_Txim_flow
  "A gas source defining mass flow, temperature and composition"
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
    powerIn=if massFlowIsLoss then 0 else min(0, gas_a.m_flow*h_port),
    powerOut=if massFlowIsLoss then 0 else max(0, gas_a.m_flow*h_port),
    powerAux=0) if                                                                                                     contributeToCycleSummary;
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean massFlowIsLoss = true
    "True if mass flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));

  parameter TILMedia.GasTypes.BaseGas                 medium = simCenter.flueGasModel
    "Medium to be used in tubes"                                                              annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  parameter Boolean variable_m_flow=false
    "True, if mass flow defined by variable input"                                       annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_T=false
    "True, if temperature defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_xi=false
    "True, if composition defined by variable input"                                      annotation(Dialog(group="Define Variable Boundaries"));

  parameter SI.MassFlowRate m_flow_const=0 "Constant mass flow rate" annotation(Dialog(group="Constant Boundaries", enable= not variable_m_flow));
  parameter SI.Temperature T_const=simCenter.T_amb_start
    "Constant specific temperature of source" annotation(Dialog(group="Constant Boundaries", enable= not variable_T));
  parameter Modelica.SIunits.MassFraction xi_const[medium.nc-1]=zeros(medium.nc-1)
    "Constant composition" annotation(Dialog(group="Constant Boundaries", enable= not variable_xi));
   TILMedia.GasObjectFunctions.GasPointer GasPointer=
        TILMedia.GasObjectFunctions.GasPointer(medium.concatGasName,8,medium.xi_default,medium.nc_propertyCalculation,medium.nc,medium.condensingIndex,0)
    "Pointer to external medium memory";

  outer ClaRa.SimCenter simCenter;
protected
  Modelica.SIunits.MassFlowRate m_flow_in;
  Modelica.SIunits.Temperature T_in;
  Modelica.SIunits.MassFraction xi_in[medium.nc-1];
  SI.EnthalpyMassSpecific h_port;

public
  ClaRa.Basics.Interfaces.GasPortOut gas_a(Medium=medium)
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  Modelica.Blocks.Interfaces.RealInput m_flow(value=m_flow_in) if (variable_m_flow)
    "Variable mass flow rate"
    annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
  Modelica.Blocks.Interfaces.RealInput T(value=T_in) if (variable_T)
    "Variable specific temperature"
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
  Modelica.Blocks.Interfaces.RealInput xi[medium.nc-1](value=xi_in) if
       (variable_xi) "Variable composition"
    annotation (Placement(transformation(extent={{-120,-80},{-80,-40}})));
equation

  h_port = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(gas_a.p,actualStream(gas_a.T_outflow),actualStream(gas_a.xi_outflow),GasPointer);

  if (not variable_m_flow) then
    m_flow_in=m_flow_const;
  end if;
  if (not variable_T) then
    T_in=T_const;
  end if;
  if (not variable_xi) then
    xi_in=xi_const;
  end if;

  gas_a.T_outflow=T_in;
  gas_a.m_flow=-m_flow_in;
  gas_a.xi_outflow=xi_in;

 annotation (Icon(graphics={
        Text(
          extent={{-100,30},{60,-30}},
          lineColor={27,36,42},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="T, xi")}),
                             Diagram(graphics));
end BoundaryGas_Txim_flow;
