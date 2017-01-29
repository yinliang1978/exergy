within Exergy.XClaRa.Components.BoundaryConditions;
model BoundaryGas_pTxi
  "A gas source defining pressure, Temperature and composition"
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

extends ClaRa.Basics.Icons.FlowSink;
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
  parameter Boolean variable_p=false
    "True, if pressure defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_T=false
    "True, if spc. temperature defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_xi=false
    "True, if composition defined by variable input"                                      annotation(Dialog(group="Define Variable Boundaries"));

  parameter SI.Pressure p_const=simCenter.p_amb_start "Constant pressure"                annotation(Dialog(group="Constant Boundaries", enable= not variable_p));
  parameter SI.Temperature T_const=simCenter.T_amb_start
    "Constant specific temperature of source"  annotation(Dialog(group="Constant Boundaries", enable= not hInputIsActive));
  parameter SI.MassFraction xi_const[medium.nc-1]=zeros(medium.nc-1)
    "Constant composition"  annotation(Dialog(group="Constant Boundaries", enable= not variable_xi));

   TILMedia.GasObjectFunctions.GasPointer GasPointer=
        TILMedia.GasObjectFunctions.GasPointer(medium.concatGasName,8,medium.xi_default,medium.nc_propertyCalculation,medium.nc,medium.condensingIndex,0)
    "Pointer to external medium memory";

  outer ClaRa.SimCenter simCenter;
protected
  SI.Pressure p_in;
  SI.Temperature T_in;
  SI.MassFraction xi_in[medium.nc-1];
  SI.EnthalpyMassSpecific h_port;

public
  Modelica.Blocks.Interfaces.RealInput p(value=p_in) if (variable_p)
    "Variable mass flow rate"
    annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
  Modelica.Blocks.Interfaces.RealInput T(value=T_in) if (variable_T)
    "Variable specific temperature"
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
  Modelica.Blocks.Interfaces.RealInput xi[medium.nc-1](value=xi_in) if
       (variable_xi) "Variable composition"
    annotation (Placement(transformation(extent={{-120,-80},{-80,-40}})));
public
  ClaRa.Basics.Interfaces.GasPortIn gas_a(Medium=medium)
    annotation (Placement(transformation(extent={{-10,80},{10,100}}),
        iconTransformation(extent={{90,-10},{110,10}})));
equation

    h_port = TILMedia.GasObjectFunctions.specificEnthalpy_pTxi(gas_a.p,actualStream(gas_a.T_outflow),actualStream(gas_a.xi_outflow),GasPointer);

  if (not variable_p) then
    p_in=p_const;
  end if;
  if (not variable_T) then
    T_in=T_const;
  end if;
  if (not variable_xi) then
    xi_in=xi_const;
  end if;

  gas_a.T_outflow=T_in;
  gas_a.p=p_in;
  gas_a.xi_outflow=xi_in;
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics={
        Text(
          extent={{-100,30},{10,-30}},
          lineColor={27,36,42},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="T
xi")}),                       Diagram(coordinateSystem(preserveAspectRatio=
            false, extent={{-100,-100},{100,100}}),
                                      graphics));
end BoundaryGas_pTxi;
