within Exergy.XClaRa.Components.HeatExchangers;
model IdealShell_L2
  "A desuperheater having an ideal cooling | block-shaped geometry"
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

  extends ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_2(redeclare model
      Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowBlockWithTubes
        (
        height=height,
        width=width,
        length=length,
        diameter_t=
            diameter_t,
        N_tubes=N_tubes,
        N_passes=N_passes,
        parallelTubes=parallelTubes,
        flowOrientation=flowOrientation,
        z_in={z_in},
        z_out={z_out}),final heatSurfaceAlloc=2, redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.IdeallyStirred);
  extends ClaRa.Basics.Icons.HEX02;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L2");
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=if not heatFlowIsLoss then -heat.Q_flow else 0,
    powerAux=0) if                                                                                                     contributeToCycleSummary;
  outer ClaRa.SimCenter   simCenter;
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // parameter dialog~~~~~~~~~~~~~~~~~
  parameter Modelica.SIunits.Length height=1 "Height of the component"
    annotation (Dialog(tab="Geometry", groupImage="modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialog_BUshell1ph2.png"));
  parameter Modelica.SIunits.Length width=1 "Width of the component"
    annotation (Dialog(tab="Geometry"));
  parameter Modelica.SIunits.Length length=1 "Length of the component"
    annotation (Dialog(tab="Geometry"));
  parameter Modelica.SIunits.Length diameter_t=0.1
    "Outer diameter of internal tubes"
    annotation (Dialog(tab="Geometry"));
  parameter Integer N_tubes=1 "Number of internal tubes"
    annotation (Dialog(tab="Geometry"));
  parameter Integer N_passes=1 "Number of passes of the internal tubes"
    annotation (Dialog(tab="Geometry"));
  parameter ClaRa.Basics.Choices.GeometryOrientation flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.horizontal
    "Flow orientation at shell side"
                                   annotation (Dialog(tab="Geometry"));
  parameter Boolean parallelTubes=false
    "True, if tubes are parallel to flow orientation, else false"
    annotation (Dialog(tab="Geometry"));

  parameter Modelica.SIunits.Length z_in=height/2 "Inlet position from bottom"
    annotation (Dialog(tab="Geometry", enable=orientation == ClaRa.Basics.Choices.GeometryOrientation.vertical));
  parameter Modelica.SIunits.Length z_out=height/2
    "Outlet position from bottom"
    annotation (Dialog(tab="Geometry", enable=orientation == ClaRa.Basics.Choices.GeometryOrientation.vertical));

  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
    annotation (Dialog(tab="Summary and Visualisation"));
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean heatFlowIsLoss = true
    "True if heat flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int
    annotation (Placement(transformation(extent={{45,-81},{47,-79}})));

public
  ClaRa.Basics.Interfaces.EyeOut eye if showData
    annotation (Placement(transformation(extent={{90,-90},{110,-70}})));

equation
  eye_int.p = outlet.p/1e5;
  eye_int.h = fluidOut.h/1e3;
  eye_int.m_flow = -outlet.m_flow;
  eye_int.T = fluidOut.T - 273.15;
  eye_int.s = fluidOut.s/1e3;
  connect(eye, eye_int) annotation (Line(
      points={{100,-80},{46,-80}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics), Diagram(graphics));
end IdealShell_L2;
