within Exergy.XClaRa.Components.BoundaryConditions;
model BoundarySlag_pT "A source defining pressure and temperature"
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
    powerIn=if massFlowIsLoss then 0 else min(0, slag_inlet.m_flow*h_slag),
    powerOut=if massFlowIsLoss then 0 else max(0, slag_inlet.m_flow*h_slag),
    powerAux=0) if                                                                                                     contributeToCycleSummary;
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean massFlowIsLoss = true
    "True if mass flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));

  parameter ClaRa.Basics.Media.Fuel.PartialSlag slagType=simCenter.slagModel
    "Medium to be used"                                                                                               annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  parameter Boolean variable_p=false
    "True, if mass flow defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_T=false
    "True, if Temperature defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
/*  parameter Boolean XInputIsActive=false 
    "True, if composition defined by variable input"                                      annotation(Dialog(group="Define Variable Boundaries"));
*/
  parameter SI.Pressure p_const=1e5 "Constant mass flow rate" annotation(Dialog(group="Constant Boundaries", enable= not mInputIsActive));
  parameter SI.Temperature T_const=simCenter.T_amb_start
    "Constant specific enthalpy of source" annotation(Dialog(group="Constant Boundaries", enable= not hInputIsActive));
 /* parameter Modelica.SIunits.MassFraction X_const[coal.nc-1]=zeros(coal.nc-1) 
    "Constant composition" annotation(Dialog(group="Constant Boundaries", enable= not XInputIsActive));
*/
  outer ClaRa.SimCenter simCenter;
protected
  Modelica.SIunits.Pressure p_in;
  Modelica.SIunits.Temperature T_in;
 // Modelica.SIunits.MassFraction X_in[coal.nc-1];
   SI.EnthalpyMassSpecific h_slag;

public
  ClaRa.Basics.Interfaces.Slag_inlet slag_inlet(final slagType=slagType)
    annotation (Placement(transformation(extent={{90,-10},{110,10}}),
        iconTransformation(extent={{90,-12},{110,8}})));

  Modelica.Blocks.Interfaces.RealInput p(value=p_in) if (variable_p)
    "Variable mass flow rate"
    annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
  Modelica.Blocks.Interfaces.RealInput T(value=T_in) if (variable_T)
    "Variable specific enthalpy"
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}})));
  /*
  Modelica.Blocks.Interfaces.RealInput X[medium.nc-1](value=X_in) if 
       (XInputIsActive) "Variable composition"
    annotation (Placement(transformation(extent={{-120,-80},{-80,-40}})));
*/
equation
  h_slag = slag_inlet.m_flow*slagType.cp*(actualStream(slag_inlet.T_outflow) - 298.15);

  if (not variable_p) then
    p_in=p_const;
  end if;
  if (not variable_T) then
    T_in=T_const;
  end if;
  /*if (not XInputIsActive) then
    X_in=X_const;
  end if;
*/
  slag_inlet.T_outflow=T_in;
  slag_inlet.p=p_in;
  //coal_a.Xi_outflow=X_in;

 annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                  graphics={
        Text(
          extent={{-100,30},{10,-30}},
          lineColor={27,36,42},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="T
xi")}),                      Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}),
                                     graphics));
end BoundarySlag_pT;
