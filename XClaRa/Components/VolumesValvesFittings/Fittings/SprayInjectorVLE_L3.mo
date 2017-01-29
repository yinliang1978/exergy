within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings;
model SprayInjectorVLE_L3 "A spray injector for i.e. temperature control"
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
  extends ClaRa.Basics.Icons.SprayInjector;
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");

  outer ClaRa.SimCenter simCenter;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=
                                                      simCenter.fluid1
    "Medium in the component" annotation(Dialog(group="Fundamental Definitions"), choicesAllMatching);
  replaceable model Material = TILMedia.SolidTypes.TILMedia_Aluminum constrainedby
    TILMedia.SolidTypes.BaseSolid "Material of the cylinder"
                               annotation(Dialog(group="Fundamental Definitions"), choicesAllMatching);

  parameter Modelica.SIunits.Length diameter_o=0.5 "Diameter of the component" annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.Length diameter_i=0.45 "Diameter of the component"
                                                                                annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.Length length=3 "Length of the component"  annotation(Dialog(group="Geometry"));
  parameter Integer N=1 "Number of identical injectors in parallel"  annotation(Dialog(group="Geometry"));

  parameter Modelica.SIunits.Pressure p_nom=1e5 "Nominal pressure" annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.SpecificEnthalpy h_nom_Main=3000e3
    "Nominal specific enthalpy"  annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.SpecificEnthalpy h_nom_Spray=1000e3
    "Nominal specific enthalpy"  annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.MassFlowRate m_flow_nom_main=300
    "Nominal mass flow rates at inlet"  annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.Pressure Delta_p_nom=1000 "Nominal pressure loss"  annotation(Dialog(group="Nominal Values"));
  parameter Modelica.SIunits.MassFlowRate m_flow_nom_spray=10
    "Nominal injection mass flow rate"  annotation(Dialog(group="Nominal Values"));

  parameter Modelica.SIunits.SpecificEnthalpy h_start_Main=3000e3
    "|Initialisation|Fluid|Initial specific enthalpy at main inlet";
  parameter Modelica.SIunits.SpecificEnthalpy h_start_Spray=1000e3
    "|Initialisation|Fluid|Initial specific enthalpy at main inlet";
  parameter Modelica.SIunits.Pressure p_start=1e5
    "|Initialisation|Fluid|Start value of sytsem pressure";

  parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "|Initialisation|Fluid|Type of initialisation";

  parameter Boolean useHomotopy=simCenter.useHomotopy
    "|Initialisation||True, if homotopy method is used during initialisation";

  ///__________Expert Settings__________________________________________________________

  parameter Boolean preciseTwoPhase = true
    "|Expert Settings|Mixing Zone|True, if two-phase transients should be calculated precisely";

  parameter Boolean checkValve = false
    "|Expert Settings|Injector Valve|True, if spray injector valve is check valve";

  parameter Boolean useStabilisedMassFlow= false
    "|Expert Settings|Injector Valve|True, if use pseudo state for mass flow in spray injector valve";

  parameter Real Tau=0.001
    "|Expert Settings|Injector Valve|Time constant of pseudo state in spray injector valve";

  parameter Real opening_leak_=0
    "|Expert Settings|Injector Valve|Leakage valve opening in p.u.";

    parameter Integer N_wall=3
    "|Expert Settings|Wall|Number of radial elements of the wall";

///__________Summary and Visualisation________________________________________________
  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "|Summary and Visualisation||True, if expert summary should be applied";
  parameter Boolean showData=false
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                            annotation(Dialog(tab="Summary and Visualisation"));
protected
  parameter Modelica.SIunits.SpecificEnthalpy h_nom_mix=(h_nom_Main*m_flow_nom_main+h_nom_Spray*m_flow_nom_spray)/(m_flow_nom_main+m_flow_nom_spray)
    "Nominal mix enthalpy";
  parameter Modelica.SIunits.SpecificEnthalpy h_start_mix=(h_start_Main*m_flow_nom_main+h_start_Spray*m_flow_nom_spray)/(m_flow_nom_main+m_flow_nom_spray)
    "Nominal mix enthalpy";
  parameter Modelica.SIunits.Temperature T_wall_start[N_wall]= ones(N_wall)*TILMedia.VLEFluidFunctions.temperature_phxi(medium, p_start-Delta_p_nom, h_start_mix)
    "Start values of wall temperature" annotation(Dialog(group="Initialisation"));

//Pressure loss models
public
  replaceable model PressureLoss =
      Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint
    constrainedby
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.GenericPressureLoss
    "Pressure loss model of injector valve"                                        annotation(Dialog(group="Fundamental Definitions"),choicesAllMatching);

  replaceable model PressureLoss_outflowZone =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L2
    "Pressure loss model of injector outflow Zone"                                                                                                     annotation(Dialog(group="Fundamental Definitions"),choicesAllMatching);

  replaceable model PressureLoss_mixingZoneInlet1 =
      Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.NoFriction
    constrainedby
    Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.BaseDp
    "Pressure loss model of mixing zone at inlet1"                                                                                        annotation(Dialog(group="Fundamental Definitions"),choicesAllMatching);

  replaceable model PressureLoss_mixingZoneInlet2 =
      Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.NoFriction
    constrainedby
    Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.BaseDp
    "Pressure loss model of mixing zone at inlet2"                                                                                        annotation(Dialog(group="Fundamental Definitions"),choicesAllMatching);

  replaceable model PressureLoss_mixingZoneOutlet =
      Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.NoFriction
    constrainedby
    Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Fundamentals.BaseDp
    "Pressure loss model of mixing zone at outlet"                                                                                        annotation(Dialog(group="Fundamental Definitions"),choicesAllMatching);

public
  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_2 outflowZone(
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.PipeGeometry (
        diameter=diameter_i,
        length=length/2),
    medium=medium,
    useHomotopy=useHomotopy,
    initType=initType,
    m_flow_nom=m_flow_nom_main + m_flow_nom_spray,
    h_nom=h_nom_mix,
    h_start=h_start_mix,
    p_nom=p_nom - Delta_p_nom,
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.IdealHeatTransfer_L2,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.IdeallyStirred,
    redeclare model PressureLoss =
        PressureLoss_outflowZone,
    p_start=p_start,
    showExpertSummary=showExpertSummary)
    annotation (Placement(transformation(extent={{40,10},{60,30}})));

  ClaRa.Basics.Interfaces.FluidPortIn inlet1(Medium=medium) "Inlet port" annotation (Placement(transformation(extent={{-110,10},{-90,30}}), iconTransformation(extent={{-110,10},{-90,30}})));
  ClaRa.Basics.Interfaces.FluidPortIn inlet2(Medium=medium) "Inlet port" annotation (Placement(transformation(extent={{-30,-110},{-10,-90}}), iconTransformation(extent={{-30,-110},{-10,-90}})));
  ClaRa.Basics.Interfaces.FluidPortOut outlet(Medium=medium) "Outlet port" annotation (Placement(transformation(extent={{90,10},{110,30}}), iconTransformation(extent={{90,10},{110,30}})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.ThickWall_L4 wall(
    redeclare replaceable model Material = Material,
    sizefunc=+1,
    diameter_o=diameter_o,
    diameter_i=diameter_i,
    length=length,
    N_tubes=N,
    initChoice=ClaRa.Basics.Choices.Init.steadyState,
    T_start=T_wall_start,
    N_rad=N_wall) annotation (Placement(transformation(extent={{40,40},{60,60}})));

  Exergy.XClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1 valve(
    medium=medium,
    openingInputIsActive=true,
    showExpertSummary=showExpertSummary,
    showData=false,
    redeclare model PressureLoss = PressureLoss,
    useStabilisedMassFlow=useStabilisedMassFlow,
    Tau=Tau,
    checkValve=checkValve,
    opening_leak_=opening_leak_) annotation (Placement(transformation(
        extent={{-10,-6},{10,6}},
        rotation=90,
        origin={-20,-46})));

  Modelica.Blocks.Interfaces.RealInput opening
    "=1: completely open, =0: completely closed" annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-40,-100}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-80,-100})));
  ClaRa.Basics.Interfaces.EyeOut eye if showData      annotation(Placement(transformation(extent={{90,-30},
            {110,-10}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int
    annotation (Placement(transformation(extent={{45,-21},{47,-19}})));
public
  Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Join_L2_Y mixingZone(
    volume=Modelica.Constants.pi/4*diameter_i^2*length/2*N,
    medium=medium,
    m_flow_in_nom={m_flow_nom_main,m_flow_nom_spray},
    p_nom=p_nom,
    h_nom=h_nom_mix,
    h_start=h_start_mix,
    p_start=p_start,
    initType=initType,
    useHomotopy=useHomotopy,
    preciseTwoPhase=preciseTwoPhase,
    showExpertSummary=showExpertSummary,
    showData=false,
    redeclare model PressureLossIn1 = PressureLoss_mixingZoneInlet1,
    redeclare model PressureLossIn2 = PressureLoss_mixingZoneInlet2,
    redeclare model PressureLossOut = PressureLoss_mixingZoneOutlet)
    annotation (Placement(transformation(extent={{-30,30},{-10,10}})));

equation
//-------------------------------------------
//Summary:
    eye_int.m_flow=-outflowZone.outlet.m_flow;
    eye_int.T= outflowZone.summary.outlet.T-273.15;
    eye_int.s=outflowZone.fluidOut.s/1e3;
    eye_int.p=outflowZone.outlet.p/1e5;
    eye_int.h=actualStream(outflowZone.outlet.h_outflow)/1e3;

  connect(wall.innerPhase, outflowZone.heat)                  annotation (Line(
      points={{49.8,40.4},{49.8,30},{50,30}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(outflowZone.outlet, outlet) annotation (Line(
      points={{60,20},{100,20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve.inlet, inlet2) annotation (Line(
      points={{-20,-56},{-20,-100}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve.opening_in, opening)     annotation (Line(
      points={{-29,-46},{-40,-46},{-40,-100}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(eye,eye_int)  annotation (Line(
      points={{100,-20},{46,-20}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(mixingZone.outlet, outflowZone.inlet) annotation (Line(
      points={{-10,20},{40,20}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(mixingZone.inlet1, inlet1) annotation (Line(
      points={{-30,20},{-100,20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(mixingZone.inlet2, valve.outlet) annotation (Line(
      points={{-20,10},{-20,-36}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,
            -100},{100,100}}),
                   graphics),    Diagram(coordinateSystem(preserveAspectRatio=true,
                  extent={{-100,-100},{100,100}})));
end SprayInjectorVLE_L3;
