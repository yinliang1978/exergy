within Exergy.XClaRa.Components.BoundaryConditions;
model BoundaryVLE_Txim_flow
  "A boundary defining temperature, composition and mass flow"
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
    powerIn=if massFlowIsLoss then 0 else min(0, steam_a.m_flow*actualStream(steam_a.h_outflow)),
    powerOut=if massFlowIsLoss then 0 else max(0, steam_a.m_flow*actualStream(steam_a.h_outflow)),
    powerAux=0) if                                                                                                     contributeToCycleSummary;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid   medium=simCenter.fluid1
    "Medium to be used"                                                                                               annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  parameter Boolean variable_m_flow=false
    "True, if mass flow defined by variable input"                                       annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_T=false
    "True, if temperature defined by variable input"                                  annotation(Dialog(group="Define Variable Boundaries"));
  parameter Boolean variable_xi=false
    "True, if composition defined by variable input"                                      annotation(Dialog(group="Define Variable Boundaries"));

  parameter SI.MassFlowRate m_flow_const=0 "Constant mass flow rate" annotation(Dialog(group="Constant Boundaries", enable= not variable_m_flow));
  parameter SI.Temperature T_const=293.15 "Constant temperature of source"
                                     annotation(Dialog(group="Constant Boundaries", enable= not variable_T));
  parameter SI.MassFraction xi_const[medium.nc-1]=zeros(medium.nc-1)
    "Constant composition" annotation(Dialog(group="Constant Boundaries", enable= not variable_xi));
  parameter SI.Pressure p_nom= 1e5 "Nominal flange pressure" annotation(Dialog(group="Nominal Values"));
  parameter SI.MassFlowRate m_flow_nom= 0
    "Nominal flange mass flow (zero refers to ideal boundary)"                                                  annotation(Dialog(group="Nominal Values"));
  outer ClaRa.SimCenter simCenter;
  parameter Boolean showData=true
    "|Summary and Visualisation||True, if a data port containing p,T,h,s,m_flow shall be shown, else false";
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean massFlowIsLoss = true
    "True if mass flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));

protected
  SI.MassFlowRate m_flow_in;
  SI.Temperature T_in;
  SI.MassFraction xi_in[medium.nc-1];
public
  ClaRa.Basics.Interfaces.FluidPortIn steam_a(Medium=medium)
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

  Modelica.Blocks.Interfaces.RealInput m_flow(value=m_flow_in) if (variable_m_flow)
    "Variable mass flow rate"
    annotation (Placement(transformation(extent={{-120,40},{-80,80}}),
        iconTransformation(extent={{-140,40},{-100,80}})));
  Modelica.Blocks.Interfaces.RealInput T(value=T_in) if (variable_T)
    "Variable temperature"
    annotation (Placement(transformation(extent={{-120,-20},{-80,20}}),
        iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealInput xi[medium.nc-1](value=xi_in) if
       (variable_xi) "Variable composition"
    annotation (Placement(transformation(extent={{-120,-80},{-80,-40}}),
        iconTransformation(extent={{-140,-80},{-100,-40}})));

  ClaRa.Basics.Interfaces.EyeOut eye if showData
    annotation (Placement(transformation(extent={{94,-86},{106,-74}})));

protected
   TILMedia.VLEFluid_pT fluidOut(
    vleFluidType=medium,
    p=steam_a.p,
    T=T_in,
    xi=xi_in) annotation (Placement(transformation(extent={{22,-20},{42,0}})));

  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{45,-81},{47,-79}}), iconTransformation(
          extent={{45,-65},{47,-63}})));

equation
  if (not variable_m_flow) then
    m_flow_in=m_flow_const;
  end if;
  if (not variable_T) then
    T_in=T_const;
  end if;
  if (not variable_xi) then
    xi_in=xi_const;
  end if;

  steam_a.h_outflow=fluidOut.h;
  if m_flow_nom>0 then
    steam_a.m_flow=-m_flow_in - (m_flow_nom/p_nom) * (p_nom - steam_a.p);
  else
    steam_a.m_flow=-m_flow_in;
  end if;
  steam_a.xi_outflow=xi_in;

  eye_int.m_flow = -steam_a.m_flow;
  eye_int.T = fluidOut.T-273.15;
  eye_int.s = fluidOut.s/1e3;
  eye_int.p = steam_a.p/1e5;
  eye_int.h = fluidOut.h/1e3;

  connect(eye,eye_int)  annotation (Line(
      points={{100,-80},{46,-80}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));

 annotation (Icon(graphics={
        Text(
          extent={{-100,30},{60,-30}},
          lineColor={27,36,42},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          textString="T, xi")}),
                             Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}),
                                     graphics));
end BoundaryVLE_Txim_flow;
