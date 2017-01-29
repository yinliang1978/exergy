within Exergy.XClaRa.Components.HeatExchangers;
model HEXvle2vle_L3_1ph_kA
  " VLE 2 VLE | L3 | 1 phase on each side | generic geometry | effective kA"
  import ClaRa;
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

  import SI = Modelica.SIunits;
  // Extends from... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  extends ClaRa.Basics.Icons.HEX01;

  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");

  outer ClaRa.SimCenter simCenter;
  model Outline
    extends ClaRa.Basics.Icons.RecordIcon;
    parameter Boolean showExpertSummary=false;
    input ClaRa.Basics.Units.HeatFlowRate Q_flow "Heat flow rate";
    input ClaRa.Basics.Units.TemperatureDifference  Delta_T_in
      "Fluid temperature at inlet T_1_in - T_2_in";
    input ClaRa.Basics.Units.TemperatureDifference  Delta_T_out
      "Fluid temperature at outlet T_1_out - T_2_out";
    input Real effectiveness if showExpertSummary "Effectivenes of HEX";
    input Real kA(unit="W/K") if showExpertSummary "Overall heat resistance";
  end Outline;

  model Summary
    extends ClaRa.Basics.Icons.RecordIcon;
    Outline outline;
  end Summary;

  //*********************************** / SHELL SIDE \ ***********************************//
  //________________________________ Shell fundamentals _______________________________//
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_shell=simCenter.fluid1
    "Medium to be used  for flow 1"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching);

  replaceable model PressureLossShell =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.ShellType_L2
    "Pressure loss model at shell side" annotation (Dialog(tab="Shell Side",
        group="Fundamental Definitions"), choicesAllMatching);

  //________________________________ Shell geometry _______________________________//
  parameter ClaRa.Basics.Units.Volume
                      volume_shell=1 "Volume of the shell side" annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length z_in_shell=1 "Inlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length z_out_shell=1
    "Outlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));

  //________________________________ Shell nominal parameter _____________________________________//

  parameter ClaRa.Basics.Units.MassFlowRate
                            m_flow_nom_shell=10
    "Nominal mass flow on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));
  parameter ClaRa.Basics.Units.Pressure
                        p_nom_shell=10 "Nominal pressure on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));
  parameter SI.SpecificEnthalpy h_nom_shell=10
    "Nominal specific enthalpy on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));

  //________________________________ Shell initialisation  _______________________________________//
  parameter SI.SpecificEnthalpy h_start_shell=1e5
    "Start value of sytsem specific enthalpy"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter ClaRa.Basics.Units.Pressure
                        p_start_shell=1e5 "Start value of sytsem pressure"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeShell=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));

  //*********************************** / TUBE SIDE \ ***********************************//
  //________________________________ Tubes fundamentals _______________________________//

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_tubes=simCenter.fluid1
    "Medium to be used for flow 2"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching);

  replaceable model PressureLossTubes =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.TubeType_L2
    "Pressure loss model at the tubes side" annotation (Dialog(tab="Tubes",
        group="Fundamental Definitions"), choicesAllMatching);

  //________________________________ Tubes geometry _______________________________//
  parameter ClaRa.Basics.Units.Volume
                      volume_tubes=1 "Volume of the tubes"     annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length z_in_tubes=1 "Inlet position from bottom"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length z_out_tubes=1
    "Outlet position from bottom"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  //________________________________ Tubes nominal parameter _____________________________________//
  parameter ClaRa.Basics.Units.MassFlowRate
                            m_flow_nom_tubes=10
    "Nominal mass flow on tube side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter ClaRa.Basics.Units.Pressure
                        p_nom_tubes=10 "Nominal pressure on tube side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter SI.SpecificEnthalpy h_nom_tubes=10
    "Nominal specific enthalpy on tube side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter ClaRa.Basics.Units.HeatFlowRate
                            Q_flow_nom=1e6 "Nominal heat flow rate"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));

  //________________________________ Tubes initialisation _______________________________________//
  parameter SI.SpecificEnthalpy h_start_tubes=1e5
    "Start value of sytsem specific enthalpy at tube side"
    annotation (Dialog(tab="Tubes", group="Initialisation"));
  parameter ClaRa.Basics.Units.Pressure
                        p_start_tubes=1e5
    "Start value of sytsem pressure at tube side"
    annotation (Dialog(tab="Tubes", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeTubes=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation at tube side"
    annotation (Dialog(tab="Tubes", group="Initialisation"));

  //*********************************** / WALL \ ***********************************//
  //________________________________ Wall fundamentals _______________________________//
  replaceable model WallMaterial =
      TILMedia.SolidTypes.TILMedia_Aluminum
    constrainedby TILMedia.SolidTypes.BaseSolid "Material of the cylinder"
    annotation (choicesAllMatching=true, Dialog(tab="Tube Wall", group=
          "Fundamental Definitions"));
  parameter ClaRa.Basics.Units.Mass  mass_struc=0
    "Mass of inner components (tubes + structural elements)"                                               annotation (Dialog(tab="Tube Wall", group="Fundamental Definitions"));
  //________________________________ Wall initialisation _______________________________________//
  parameter ClaRa.Basics.Units.Temperature T_w_i_start=293.15
    "Initial temperature at inner phase"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));
  parameter ClaRa.Basics.Units.Temperature T_w_o_start=293.15
    "Initial temperature at outer phase"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initWall=ClaRa.Basics.Choices.Init.noInit
    "Initialisation option"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));

//*******************************General *************************************************//
  parameter Boolean tubesLimitHeatFlow = true
    "True if the tube side heat transfer limits overall performance"                                           annotation(Dialog(tab="General", group="Heat Exchanger Definition"));
  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching);
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  parameter Real kA(unit="W/K")=50000 "The product kA - nominal value" annotation(Dialog(tab="General", group="Heat Exchanger Definition"));
  parameter Real CL_kA_mflow[:, 2]=[0,0;0.4, 0.5; 0.5, 0.75; 0.75, 0.95;
      1, 1] "Characteristic line kA = f(m_flow/m_flow_nom)" annotation(Dialog(tab="General", group="Heat Exchanger Definition"));
  replaceable model HeatExchangerType =
      ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.HeatExchangerTypes.CounterFlow
    constrainedby
    ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.HeatExchangerTypes.GeneralHeatExchanger
    "Heat exchanger type"                                                                     annotation(choicesAllMatching, Dialog(tab="General", group="Heat Exchanger Definition"));
//   parameter Real CF_geo=1
//     "|Heat Exchanger Definition|Correction coefficient due to fins etc.";

  //*********************************** / EXPERT Settings and Visualisation \ ***********************************//

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "|Summary and Visualisation||True, if expert summary should be applied";
  parameter Boolean showData=true
    "|Summary and Visualisation||True, if a data port containing p,T,h,s,m_flow shall be shown, else false";

  ClaRa.Basics.Interfaces.FluidPortIn In2(Medium=medium_tubes)
    annotation (Placement(transformation(extent={{90,-70},{110,-50}}), iconTransformation(extent={{90,-70},{110,-50}})));
  ClaRa.Basics.Interfaces.FluidPortOut Out2(Medium=medium_tubes)
    annotation (Placement(transformation(extent={{90,50},{110,70}}),
        iconTransformation(extent={{90,50},{110,70}})));
  ClaRa.Basics.Interfaces.FluidPortOut Out1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
  ClaRa.Basics.Interfaces.FluidPortIn In1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,88},{10,108}})));

  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_2 tubes(
    final heatSurfaceAlloc=1,
    medium=medium_tubes,
    p_nom=p_nom_tubes,
    h_nom=h_nom_tubes,
    m_flow_nom=m_flow_nom_tubes,
    useHomotopy=useHomotopy,
    h_start=h_start_tubes,
    p_start=p_start_tubes,
    initType=initTypeTubes,
    redeclare model PressureLoss = PressureLossTubes,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.IdeallyStirred,
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.IdealHeatTransfer_L2,
    showExpertSummary=showExpertSummary,
    redeclare model Geometry =
        Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry (
        volume=volume_tubes,
        z_in={z_in_tubes},
        z_out={z_out_tubes},
        N_heat=1))
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={70,0})));

  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_2 shell(
    medium=medium_shell,
    p_nom=p_nom_shell,
    h_nom=h_nom_shell,
    redeclare model PressureLoss = PressureLossShell,
    m_flow_nom=m_flow_nom_shell,
    useHomotopy=useHomotopy,
    h_start=h_start_shell,
    p_start=p_start_shell,
    initType=initTypeShell,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.IdeallyStirred,
    showExpertSummary=showExpertSummary,
    final heatSurfaceAlloc=1,
    redeclare model HeatTransfer =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.IdealHeatTransfer_L2,
    redeclare model Geometry =
        Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry (
        volume=volume_shell,
        z_in={z_in_shell},
        z_out={z_out_shell},
        N_heat=1))    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,0})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.NTU_L2_effectiveResistance wall(
    mass_struc=mass_struc,
    T_i_in=tubes.fluidIn.T,
    T_a_in=shell.fluidIn.T,
    m_flow_i=tubes.inlet.m_flow,
    m_flow_a=shell.inlet.m_flow,
    cp_mean_i=(tubes.fluidIn.cp + tubes.fluidOut.cp)/2,
    cp_mean_a=(shell.fluidIn.cp + shell.fluidOut.cp)/2,
    kA=kA,
    redeclare model Material = WallMaterial,
    CL_kA_mflow=CL_kA_mflow,
    redeclare model HeatExchangerType = HeatExchangerType,
    T_w_i_start=T_w_i_start,
    T_w_a_start=T_w_o_start,
    initChoice=initWall,
    m_flow_nom=if tubesLimitHeatFlow then m_flow_nom_tubes else
        m_flow_nom_shell,
    innerSideLimitsHeatFlow=tubesLimitHeatFlow) annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={30,0})));

  Summary summary(outline(
      showExpertSummary=showExpertSummary,
      Q_flow=shell.heat.Q_flow,
      Delta_T_in=shell.summary.inlet.T - tubes.summary.inlet.T,
      Delta_T_out=shell.summary.outlet.T - tubes.summary.outlet.T,
      effectiveness=wall.effectiveness,
      kA=wall.kA*wall.partLoad_kA.y[1])) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-50,-92})));
  ClaRa.Basics.Interfaces.EyeOut eye2 if showData annotation (Placement(
        transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={100,-86}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={110,-90})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int2
    annotation (Placement(transformation(extent={{63,-87},{65,-85}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int1
    annotation (Placement(transformation(extent={{27,-59},{29,-57}})));
public
  ClaRa.Basics.Interfaces.EyeOut eye1 if showData annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={30,-110}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={30,-110})));
initial equation
  //        wall.T=(tubes.bulk.T+shell.bulk.T)/2;

equation
  eye_int1.m_flow = shell.summary.outlet.m_flow;
  eye_int1.T=shell.summary.outlet.T-273.15;
  eye_int1.s=shell.fluidOut.s/1000;
  eye_int1.h=shell.summary.outlet.h/1000;
  eye_int1.p=shell.summary.outlet.p/100000;

eye_int2.m_flow = tubes.summary.outlet.m_flow;
  eye_int2.T=tubes.summary.outlet.T-273.15;
  eye_int2.s=tubes.fluidOut.s/1000;
  eye_int2.h=tubes.summary.outlet.h/1000;
  eye_int2.p=tubes.summary.outlet.p/100000;
  connect(tubes.inlet, In2) annotation (Line(
      points={{70,-10},{70,-60},{100,-60}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(tubes.outlet, Out2) annotation (Line(
      points={{70,10},{70,60},{100,60}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shell.inlet, In1) annotation (Line(
      points={{1.77636e-015,10},{1.77636e-015,74},{0,74},{0,98}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shell.outlet, Out1) annotation (Line(
      points={{-1.77636e-015,-10},{0,-10},{0,-100}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(wall.innerPhase, tubes.heat) annotation (Line(
      points={{39,-4.44089e-016},{39,4.44089e-016},{60,4.44089e-016}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(wall.outerPhase, shell.heat) annotation (Line(
      points={{21,4.44089e-016},{21,-1.77636e-015},{10,-1.77636e-015}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(eye_int1, eye1) annotation (Line(points={{28,-58},{30,-58},{30,-110}}, color={190,190,190}));
  connect(eye2, eye_int2) annotation (Line(points={{100,-86},{100,-86},{64,-86}},
                                                                                color={190,190,190}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics={Text(
          extent={{-90,94},{82,54}},
          lineColor={27,36,42},
          textString="NTU")}),Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end HEXvle2vle_L3_1ph_kA;
