within Exergy.XClaRa.Components.MechanicalSeparation;
partial model FeedWaterTank_base "Base class for feedwater tanks"
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                            //
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
  extends ClaRa.Basics.Icons.FeedwaterTank;

  outer ClaRa.SimCenter simCenter;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid
                                      medium=simCenter.fluid1
    "Medium in the component" annotation(Dialog(group="Fundamental Definitions"), choicesAllMatching);

  parameter Modelica.SIunits.Length diameter=1 "Diameter of the component"  annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.Length length=1 "Length of the component"  annotation(Dialog(group="Geometry"));
  parameter ClaRa.Basics.Choices.GeometryOrientation orientation=ClaRa.Basics.Choices.GeometryOrientation.vertical
    "Orientation of the component" annotation (Dialog(group="Geometry"));

  parameter Modelica.SIunits.MassFlowRate m_flow_cond_nom
    "Nominal condensat flow"                                                       annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.MassFlowRate m_flow_heat_nom
    "Nominal heating steam flow"                                                       annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.Pressure p_nom=1e5 "Nominal pressure"  annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.SpecificEnthalpy h_nom=1e5
    "Nominal specific enthalpy"  annotation(Dialog(group="Nominal Values"));
  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation"));

  parameter Modelica.SIunits.Pressure p_start=1e5
    "Start value of sytsem pressure"                                                   annotation(Dialog(tab="Initialisation"));
  parameter Real level_rel_start "Initial filling level" annotation(Dialog(tab="Initialisation"));
  final parameter ClaRa.Basics.Units.MassFraction xi_start=(1 - level_rel_start)*
      TILMedia.VLEFluidFunctions.dewDensity_pxi(medium, p_start)/((1 - level_rel_start)
      *TILMedia.VLEFluidFunctions.dewDensity_pxi(medium, p_start) + level_rel_start*
      TILMedia.VLEFluidFunctions.bubbleDensity_pxi(medium, p_start));
  parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation" annotation (Dialog(tab="Initialisation"));

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "True, if expert summary should be applied"                                                               annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown"                              annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean showLevel = false "True, if level shall be visualised"  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean levelOutput = false
    "True, if Real level connector shall be addded"                                      annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean outputAbs = false "True, if absolute level is at output"  annotation(Dialog(enable = levelOutput, tab="Summary and Visualisation"));

public
  ClaRa.Basics.Interfaces.FluidPortOut outlet(Medium=medium) "Outlet port"
    annotation (Placement(transformation(extent={{-270,-110},{-250,-90}}), iconTransformation(extent={{-270,-110},{-250,-90}})));
  ClaRa.Basics.Interfaces.FluidPortIn heatingSteam(Medium=medium)
    "Heating steam inlet"
    annotation (Placement(transformation(extent={{-210,70},{-190,90}}), iconTransformation(extent={{-210,70},{-190,90}})));
  ClaRa.Basics.Interfaces.FluidPortIn condensate(Medium=medium)
    "Main condensate inlet"
    annotation (Placement(transformation(extent={{190,50},{210,70}}), iconTransformation(extent={{190,50},{210,70}})));
public
  ClaRa.Basics.Interfaces.EyeOut eye if showData
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-220,-110}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-220,-110})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int
    annotation (Placement(transformation(extent={{-221,-81},{-219,-79}})));
equation
  connect(eye_int, eye) annotation (Line(
      points={{-220,-80},{-220,-110}},
      color={190,190,190},
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-300,-100},{300,100}},
        initialScale=0.1)),
                         Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-300,-100},{300,100}},
        initialScale=0.1)));
end FeedWaterTank_base;
