within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings;
model SprayInjectorVLE_L3_advanced
  "A spray injector for i.e. temperature control"
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

///_______________Fundamental Definitions__________________________________________________
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium=
                                                      simCenter.fluid1
    "Medium in the component" annotation(Dialog(group="Fundamental Definitions"), choicesAllMatching);

  replaceable model Material = TILMedia.SolidTypes.TILMedia_Aluminum constrainedby
    TILMedia.SolidTypes.BaseSolid "Material of the cylinder"
                               annotation(Dialog(group="Fundamental Definitions"), choicesAllMatching);

  replaceable model PressureLoss =
      Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint
    constrainedby
    Exergy.XClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.GenericPressureLoss
    "Pressure loss model of injector valve"                                        annotation(Dialog(group="Fundamental Definitions"),choicesAllMatching);

  replaceable model PressureLoss_mixingZone =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.QuadraticParallelZones_L3
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L3
    "Pressure loss model of injector mixing Zone"                                                                                                     annotation(Dialog(group="Fundamental Definitions"),choicesAllMatching);

  parameter Modelica.SIunits.Length diameter_o=0.5 "Diameter of the component" annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.Length diameter_i=0.45 "Diameter of the component"
                                                                                annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.Length length=3 "Length of the component"  annotation(Dialog(group="Geometry"));

  parameter Modelica.SIunits.Pressure p_nom=1e5 "Nominal pressure" annotation(Dialog(group="Nominal Values"));
//     parameter Modelica.SIunits.SpecificEnthalpy h_nom_Main=3000e3 "Nominal specific enthalpy"
//                                    annotation(Dialog(group="Nominal Values"));
//    parameter Modelica.SIunits.SpecificEnthalpy h_nom_Spray=1000e3 "Nominal specific enthalpy"
//                                   annotation(Dialog(group="Nominal Values"));
    parameter Modelica.SIunits.MassFlowRate m_flow_nom_main=10
    "Nominal mass flow rates at inlet"    annotation(Dialog(group="Nominal Values"));
//
//    parameter Modelica.SIunits.MassFlowRate m_flow_nomSpray=10 "Nominal injection mass flow rate"
//                                          annotation(Dialog(group="Nominal Values"));

///__________INitialisation____________________________________________________________
   parameter Modelica.SIunits.SpecificEnthalpy h_start_main=3000e3
    "|Initialisation|Fluid|Initial specific enthalpy of main phase";
   parameter Modelica.SIunits.SpecificEnthalpy h_start_spray=h_start_main-100
    "|Initialisation|Fluid|Initial specific enthalpy of spray phase";
  parameter Modelica.SIunits.Pressure p_start=1e5
    "|Initialisation|Fluid|Start value of sytsem pressure";

  parameter Real y_start=0.05
    "|Initialisation|Fluid|Start value for ratio spray volume to total volume";
  parameter ClaRa.Basics.Choices.Init initFluid=ClaRa.Basics.Choices.Init.noInit
    "|Initialisation|Fluid|Initialisation option of fluid"
                                              annotation(Dialog(group="Initialisation"));

   parameter Modelica.SIunits.Temperature T_wall_start[N_wall]=ones(N_wall)*TILMedia.VLEFluidFunctions.temperature_phxi(medium, p_start, h_start_main)
    "|Initialisation|Wall|Start values of wall temperature";
  parameter ClaRa.Basics.Choices.Init initWall=ClaRa.Basics.Choices.Init.noInit
    "|Initialisation|Wall|Initialisation option of wall";

  parameter Boolean useHomotopy=simCenter.useHomotopy
    "|Initialisation||True, if homotopy method is used during initialisation";

///__________Expert Settings__________________________________________________________
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

  parameter SI.Time Tau_cond=0.03
    "|Expert Settings|Mixing Model|Time constant of condensation";
  parameter SI.Time Tau_evap=Tau_cond
    "|Expert Settings|Mixing Model|Time constant of evaporation";
  parameter SI.CoefficientOfHeatTransfer alpha_ph=50000
    "|Expert Settings|Mixing Model|HTC of the phase border";
  parameter SI.Area A_phaseBorder=10
    "|Expert Settings|Mixing Model|Heat transfer area at phase border";
  parameter Real exp_HT_phases=0
    "|Expert Settings|Mixing Model|Exponent for volume dependency on inter phase HT";

  parameter SI.VolumeFraction eps_mix[2]={0.2,0.8}
    "|Expert Settings|Mixing Model|Volume fraction V_1/V_tot of min/max mixed outlet";

///__________Summary and Visualisation________________________________________________
  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "|Summary and Visualisation||True, if expert summary should be applied";

  parameter Boolean showData=false
    "|Summary and Visualisation||True, if a data port containing p,T,h,s,m_flow shall be shown, else false";

//   parameter Modelica.SIunits.SpecificEnthalpy h_nom_mix=(h_nom_Main*m_flow_nom_main+h_nom_Spray*m_flow_nomSpray)/(m_flow_nom_main+m_flow_nomSpray)
//     "Nominal mix enthalpy";

public
  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_L3_TwoZonesNPort
    mixingZone(
    medium=medium,
    useHomotopy=useHomotopy,
    initType=initFluid,
    p_start=p_start,
    p_nom=p_nom,
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry (
        N_inlet=2,
        A_cross=Modelica.Constants.pi/4*diameter_i^2,
        A_front=Modelica.Constants.pi/4*diameter_i^2,
        A_hor=diameter_i*length,
        volume=Modelica.Constants.pi/4*diameter_i^2*length,
        N_heat=1,
        A_heat={Modelica.Constants.pi*diameter_i*length},
        height_fill=diameter_i),
    m_flow_nom=m_flow_nom_main,
    Tau_cond=Tau_cond,
    Tau_evap=Tau_evap,
    alpha_ph=alpha_ph,
    A_heat_ph=A_phaseBorder,
    exp_HT_phases=exp_HT_phases,
    h_liq_start=h_start_spray,
    h_vap_start=h_start_main,
    level_rel_start=y_start,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.RealMixed
        (level_rel_start=y_start, eps_mix=eps_mix),
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3,
    showExpertSummary=showExpertSummary,
    redeclare model PressureLoss = PressureLoss_mixingZone)
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
    T_start=T_wall_start,
    N_rad=N_wall,
    N_tubes=1,
    initChoice=initWall) annotation (Placement(transformation(extent={{0,60},{20,80}})));

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

  Adapters.Scalar2VectorHeatPort scalar2VectorHeatPort(N=2)
    annotation (Placement(transformation(extent={{20,40},{40,60}})));
equation
//-------------------------------------------
//Summary:
    eye_int.m_flow=-mixingZone.outlet[1].m_flow;
    eye_int.T= mixingZone.summary.outlet[1].T-273.15;
    eye_int.s=mixingZone.fluidOut[1].s/1e3;
    eye_int.p=mixingZone.outlet[1].p/1e5;
    eye_int.h=actualStream(mixingZone.outlet[1].h_outflow)/1e3;

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
  connect(inlet1, mixingZone.inlet[1]) annotation (Line(
      points={{-100,20},{40,20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valve.outlet, mixingZone.inlet[2])  annotation (Line(
      points={{-20,-36},{-20,20},{40,20}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(outlet, mixingZone.outlet[1]) annotation (Line(
      points={{100,20},{60,20}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(scalar2VectorHeatPort.heatVector, mixingZone.heat) annotation (Line(
      points={{40,50},{50,50},{50,30}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(scalar2VectorHeatPort.heatScalar, wall.innerPhase) annotation (Line(
      points={{20,50},{10,50},{10,60.4},{9.8,60.4}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},
            {100,100}}),
                   graphics),    Diagram(coordinateSystem(preserveAspectRatio=true,
                  extent={{-100,-100},{100,100}})));
end SprayInjectorVLE_L3_advanced;
