within Exergy.XClaRa.Components.HeatExchangers;
model HEXvle2vle_L3_2ph_CH_simple
  "VLE 2 VLE | L3 | 2 phase at shell side | Cylinder shape | Header type | simple HT"
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

  // Extends from... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  extends ClaRa.Basics.Icons.HEX05;

  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");

  outer ClaRa.SimCenter simCenter;

  model Outline
    extends ClaRa.Basics.Icons.RecordIcon;
    parameter Boolean showExpertSummary=false;
    input SI.HeatFlowRate Q_flow "Heat flow rate";
    input SI.TemperatureDifference Delta_T_in
      "Fluid temperature at inlet T_1_in - T_2_in";
    input SI.TemperatureDifference Delta_T_out
      "Fluid temperature at outlet T_1_out - T_2_out";
    input SI.Length absLevel "Absolute filling level";
    input Real relLevel "relative filling level";
  end Outline;

  model Summary
    extends ClaRa.Basics.Icons.RecordIcon;
    Outline outline;
  end Summary;

  // Parameters and other user definable settings~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  //*********************************** / SHELL SIDE \ ***********************************//
  //________________________________ Shell fundamentals _______________________________//
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_shell=simCenter.fluid1
    "Medium to be used for shell flow" annotation (choices(
      choice=simCenter.fluid1 "First fluid defined in global simCenter",
      choice=simCenter.fluid2 "Second fluid defined in global simCenter",
      choice=simCenter.fluid3 "Third fluid defined in global simCenter"),
                                                          Dialog(tab=
          "Shell Side", group="Fundamental Definitions"));
  replaceable model HeatTransfer_Shell =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.ShellType_L3
    "Heat transfer model at shell side" annotation (Dialog(tab="Shell Side",
        group="Fundamental Definitions"), choicesAllMatching);
  replaceable model PressureLossShell =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L3
    "Pressure loss model at shell side" annotation (Dialog(tab="Shell Side",
        group="Fundamental Definitions"), choicesAllMatching);
  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"
   annotation (Dialog(tab="General",group="Fundamental Definitions"), choicesAllMatching);
//, groupImage="modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialog_CHgeneral.png"
  //________________________________ Shell geometry _______________________________//
  parameter SI.Length length=10 "Length of the HEX" annotation (Dialog(
      tab="Shell Side",
      group="Geometry",
      groupImage=
          "modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialog_CHshell.png"));                                                             //
  parameter SI.Length diameter=3 "Diameter of HEX"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter SI.Length z_in_shell=length/2 "Inlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
    parameter SI.Length z_in_aux1=length/2
    "Inlet position of auxilliary1 from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter SI.Length z_in_aux2=length/2
    "Inlet position of auxilliary2 from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter SI.Length z_out_shell=length/2 "Outlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter SI.Length radius_flange=0.05 "Flange radius of all flanges"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  final parameter SI.Mass mass_struc=0
    "Mass of inner structure elements, additional to the tubes itself"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Choices.GeometryOrientation orientation=ClaRa.Basics.Choices.GeometryOrientation.vertical
    "Orientation of the component"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Choices.GeometryOrientation flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.vertical
    "Orientation of the mass flow"
    annotation (Dialog(tab="Shell Side", group="Geometry"));

  //________________________________ Shell nominal parameter _____________________________________//
  parameter SI.MassFlowRate m_flow_nom_shell=10
    "Nominal mass flow on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));
  parameter SI.Pressure p_nom_shell=10 "Nominal pressure on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));
  parameter SI.EnthalpyMassSpecific h_nom_shell=100e3
    "Nominal specific enthalpy on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));

  //________________________________ Shell initialisation  _______________________________________//
    parameter SI.EnthalpyMassSpecific h_liq_start=-10 +
      TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(medium_shell,
      p_start_shell) "Start value of liquid specific enthalpy"
                                                              annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter SI.EnthalpyMassSpecific h_vap_start=+10 +
      TILMedia.VLEFluidFunctions.dewSpecificEnthalpy_pxi(medium_shell, p_start_shell)
    "Start value of vapour specific enthalpy"                                                                                   annotation (Dialog(tab="Shell Side", group="Initialisation"));

  parameter SI.Pressure p_start_shell=1e5 "Start value of shell fluid pressure"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter Real level_rel_start=0.5 "Start value for relative filling Level" annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeShell=ClaRa.Basics.Choices.Init.noInit
    "Type of shell fluid initialisation"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));

  //*********************************** / TUBE SIDE \ ***********************************//
  //________________________________ Tubes fundamentals _______________________________//

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_tubes=simCenter.fluid1
    "Medium to be used for tubes flow"                                                                           annotation (choices(
      choice=simCenter.fluid1 "First fluid defined in global simCenter",
      choice=simCenter.fluid2 "Second fluid defined in global simCenter",
      choice=simCenter.fluid3 "Third fluid defined in global simCenter"),
                                                          Dialog(tab="Tubes",
        group="Fundamental Definitions"));
  replaceable model HeatTransferTubes =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.CharLine_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.TubeType_L2
    "Heat transfer mode at the tubes side"                                                                                  annotation (Dialog(tab="Tubes",
        group="Fundamental Definitions"), choicesAllMatching);
  replaceable model PressureLossTubes =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.TubeType_L2
    "Pressure loss model at the tubes side" annotation (Dialog(tab="Tubes",group="Fundamental Definitions"), choicesAllMatching);

  //________________________________ Tubes geometry _______________________________//
  parameter SI.Length diameter_i=0.048 "Inner diameter of internal tubes"
                                       annotation (Dialog(
      tab="Tubes",
      group="Tubes Geometry",
      groupImage=
          "modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialogTubes.png"));
  parameter SI.Length diameter_o=0.05 "Outer diameter of internal tubes"
    annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter SI.Length length_tubes=10 "Length of the tubes (one pass)"
    annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter Integer N_tubes=1000 "Number of tubes"  annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter Integer N_passes=1 "Number of passes of the internal tubes" annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter Boolean parallelTubes=true
    "True, if tubes are parallel to shell flow orientation"                                    annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter SI.Length z_in_tubes=length/2 "Inlet position from bottom"
    annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter SI.Length z_out_tubes=length/2 "Outlet position from bottom"
    annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter Boolean staggeredAlignment=true
    "True, if the tubes are aligned staggeredly"                                         annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter SI.Length Delta_z_par=2*diameter_o
    "Distance between tubes parallel to flow direction (center to center)"
    annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter SI.Length Delta_z_ort=2*diameter_o
    "Distance between tubes orthogonal to flow direction (center to center)"
    annotation (Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter Integer N_rows=integer(ceil(sqrt(N_tubes))*N_passes)
    "Number of pipe rows in shell flow direction"                                                              annotation(Dialog(tab="Tubes", group="Tubes Geometry"));
  parameter Real CF_geo=1 "Correction coefficient due to fins etc." annotation (Dialog(tab="Tubes", group="Tubes Geometry"));

  //________________________________ Tubes nominal parameter _____________________________________//
  parameter SI.MassFlowRate m_flow_nom_tubes=10
    "Nominal mass flow on tubes side" annotation (Dialog(
      tab="Tubes",
      group="Nominal Values",
      groupImage=
          "modelica://ClaRa/figures/ParameterDialog/CH_general.png"));
  parameter SI.Pressure p_nom_tubes=10 "Nominal pressure on side tubes"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter SI.EnthalpyMassSpecific h_nom_tubes=10
    "Nominal specific enthalpy on tubes side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter SI.HeatFlowRate Q_flow_nom=1e6 "Nominal heat flow rate"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));

  //___________________________Initialisation tubes _______________________________________//
  parameter SI.EnthalpyMassSpecific h_start_tubes=1e5
    "Start value of tube fluid specific enthalpy"
    annotation (Dialog(tab="Tubes", group="Initialisation"));
  parameter SI.Pressure p_start_tubes=1e5 "Start value of tube fluid pressure"
    annotation (Dialog(tab="Tubes", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeTubes=ClaRa.Basics.Choices.Init.noInit
    "Type of tube fluid initialisation"
    annotation (Dialog(tab="Tubes", group="Initialisation"));

  //***********************************/ WALL \ *****************************************//
  replaceable model WallMaterial =
      TILMedia.SolidTypes.TILMedia_Aluminum
    constrainedby TILMedia.SolidTypes.BaseSolid "Material of the cylinder"
    annotation (choicesAllMatching=true, Dialog(tab="Tube Wall", group=
          "Fundamental Definitions"));

  parameter SI.Temperature T_w_start[3]=ones(3)*293.15
    "Initial wall temperature inner --> outer"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeWall=ClaRa.Basics.Choices.Init.noInit
    "Init option of Tube wall"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));

  //*********************************** / EXPERT Settings and Visualisation \ ***********************************//
  parameter SI.Time Tau_cond=0.3 "Time constant of condensation"
    annotation (Dialog(tab="Expert Settings", group=
          "Zone Interaction at Shell Side"));
  parameter SI.Time Tau_evap=0.03 "Time constant of evaporation"
    annotation (Dialog(tab="Expert Settings", group=
          "Zone Interaction at Shell Side"));
  parameter SI.CoefficientOfHeatTransfer alpha_ph=50000
    "HTC of the phase border" annotation (Dialog(tab="Expert Settings",
        group="Zone Interaction at Shell Side"));
  parameter SI.Area A_phaseBorder=shell.geo.A_hor*100
    "Heat transfer area at phase border" annotation (Dialog(tab=
          "Expert Settings", group="Zone Interaction at Shell Side"));
  parameter Real expHT_phases=0
    "Exponent for volume dependency on inter phase HT"                             annotation (Dialog(tab="Expert Settings", group="Zone Interaction at Shell Side"));
  parameter Real absorbInflow=1
    "Absorption of incoming mass flow to the zones 1: perfect in the allocated zone, 0: perfect according to steam quality"
                                                                                              annotation (Dialog(tab="Expert Settings", group="Zone Interaction at Shell Side"));
  parameter Boolean equalPressures=true
    "True if pressure in liquid and vapour phase is equal"                                     annotation (Dialog(tab="Expert Settings", group="Zone Interaction at Shell Side"));
 parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
    "Smoothness of level calculation (table based)"                                                                                                    annotation (Dialog(tab="Expert Settings", group="Mass Accumulation at Shell Side"));

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "True, if expert summary should be applied"                                                               annotation (Dialog(tab="Summary and Visualisation"));
  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                              annotation (Dialog(tab="Summary and Visualisation"));
  parameter Boolean levelOutput = false
    "True, if Real level connector shall be addded"                                      annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean outputAbs = false "True, if absolute level is at output"  annotation(Dialog(enable = levelOutput, tab="Summary and Visualisation"));

  ClaRa.Basics.Interfaces.FluidPortIn In2(Medium=medium_tubes)
    annotation (Placement(transformation(extent={{90,10},{110,30}}),
        iconTransformation(extent={{90,10},{110,30}})));
  ClaRa.Basics.Interfaces.FluidPortOut Out2(Medium=medium_tubes)
    annotation (Placement(transformation(extent={{-110,10},{-90,30}}),
        iconTransformation(extent={{-110,10},{-90,30}})));
  ClaRa.Basics.Interfaces.FluidPortOut Out1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
  ClaRa.Basics.Interfaces.FluidPortIn In1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,88},{10,108}})));

  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_2 tubes(
    medium=medium_tubes,
    p_nom=p_nom_tubes,
    h_nom=h_nom_tubes,
    m_flow_nom=m_flow_nom_tubes,
    useHomotopy=useHomotopy,
    h_start=h_start_tubes,
    p_start=p_start_tubes,
    initType=initTypeTubes,
    redeclare model HeatTransfer = HeatTransferTubes,
    redeclare model PressureLoss = PressureLossTubes,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.IdeallyStirred,
    showExpertSummary=showExpertSummary,
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.PipeGeometry (
        final z_in={z_in_tubes},
        final z_out={z_out_tubes},
        final diameter=diameter_i,
        final N_tubes=N_tubes,
        final N_passes=N_passes,
        final length=length_tubes))
    annotation (Placement(transformation(extent={{62,10},{42,30}})));

  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_L3_TwoZonesNPort shell(
    medium=medium_shell,
    p_nom=p_nom_shell,
    redeclare model HeatTransfer = HeatTransfer_Shell,
    redeclare model PressureLoss = PressureLossShell,
    m_flow_nom=m_flow_nom_shell,
    useHomotopy=useHomotopy,
    p_start=p_start_shell,
    initType=initTypeShell,
    level_rel_start=level_rel_start,
    Tau_cond=Tau_cond,
    Tau_evap=Tau_evap,
    alpha_ph=alpha_ph,
    h_liq_start=h_liq_start,
    h_vap_start=h_vap_start,
    showExpertSummary=showExpertSummary,
    exp_HT_phases=expHT_phases,
    A_heat_ph=A_phaseBorder,
    heatSurfaceAlloc=2,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.RealSeparated
        (
        level_rel_start=level_rel_start,
        radius_flange=radius_flange,
        absorbInflow=absorbInflow,
        smoothness=smoothness),
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowCylinderWithTubes
        (
        final z_out={z_out_shell},
        final length=length,
        final N_tubes=N_tubes,
        final staggeredAlignment=staggeredAlignment,
        final Delta_z_par=Delta_z_par,
        final Delta_z_ort=Delta_z_ort,
        final diameter=diameter,
        final N_passes=N_passes,
        final length_tubes=length_tubes,
        final diameter_t=diameter_o,
        final N_inlet=3,
        final z_in={z_in_shell,z_in_aux1,z_in_aux2},
        final orientation=orientation,
        final flowOrientation=flowOrientation,
        final N_rows=N_rows,
        final parallelTubes=parallelTubes,
        final N_baffle=0,
        CF_geo={1,CF_geo}),
    equalPressures=equalPressures)
                           annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,46})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.ThickWall_L4 wall(
    redeclare replaceable model Material = WallMaterial,
    initChoice=initTypeWall,
    N_rad=3,
    sizefunc=1,
    diameter_o=diameter_o,
    diameter_i=diameter_i,
    N_tubes=N_tubes,
    T_start=T_w_start,
    length=length*N_passes) annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={33,45})));

protected
TILMedia.VLEFluid_ph fluidIn1(
    each vleFluidType=medium_shell,
    h=actualStream(In1.h_outflow),
    p=In1.p)                                                           annotation (Placement(transformation(extent={{10,72},
            {30,92}},                                                                                                   rotation=0)));
TILMedia.VLEFluid_ph fluidOut1(
    each vleFluidType=medium_shell,
    h=actualStream(Out1.h_outflow),
    p=Out1.p)                                                         annotation (Placement(transformation(extent={{2,-70},
            {22,-50}},                                                                                                  rotation=0)));

 TILMedia.VLEFluid_ph fluidIn2(
    each vleFluidType=medium_tubes,
    h=actualStream(In2.h_outflow),
    p=In2.p)                                                           annotation (Placement(transformation(extent={{72,-12},
            {92,8}},                                                                                                    rotation=0)));
TILMedia.VLEFluid_ph fluidOut2(
    each vleFluidType=medium_tubes,
    h=actualStream(Out2.h_outflow),
    p=Out2.p)                                                         annotation (Placement(transformation(extent={{-90,-12},
            {-70,8}},                                                                              rotation=0)));

public
    outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Exergy.Utilities.ViewObjectNE viewObject(nEnergy={4,1,0,0});

  Exergy.Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));
equation

  //  viewObject.h[1].E_flow = inlet.m_flow*noEvent(actualStream(inlet.h_outflow));
  viewObject.h[1].E_flow = In1.m_flow*fluidIn1.h;
  viewObject.h[1].Ex_flow = In1.m_flow*(fluidIn1.h-refEnv.T*fluidIn1.s);

  viewObject.h[2].E_flow = Out1.m_flow*fluidOut1.h;
  viewObject.h[2].Ex_flow = Out1.m_flow*(fluidOut1.h-refEnv.T*fluidOut1.s);

    //  viewObject.h[1].E_flow = inlet.m_flow*noEvent(actualStream(inlet.h_outflow));
  viewObject.h[3].E_flow = In2.m_flow*fluidIn2.h;
  viewObject.h[3].Ex_flow = In2.m_flow*(fluidIn2.h-refEnv.T*fluidIn2.s);

  viewObject.h[4].E_flow = Out2.m_flow*fluidOut2.h;
  viewObject.h[4].Ex_flow = Out2.m_flow*(fluidOut2.h-refEnv.T*fluidOut2.s);

    viewObject.E[1].E = shell.summary.fluid.mass * shell.summary.fluid.h
                       +tubes.mass * tubes.bulk.h;
   viewObject.E[1].Ex = shell.summary.fluid.mass * (shell.summary.fluid.h - refEnv.T*{shell.liq.s,shell.vap.s})
                          +tubes.mass * (tubes.bulk.h - refEnv.T*tubes.bulk.s);
   //viewObject.E[1].Ex = shell.summary.fluid.mass * ( refEnv.T*(shell.summary.fluid.rho));

  connect(viewObject.viewOutput,viewOutput);

public
  Summary summary(outline(
      showExpertSummary=showExpertSummary,
      Q_flow=sum(shell.heat.Q_flow),
      Delta_T_in=shell.summary.inlet[1].T - tubes.summary.inlet.T,
      Delta_T_out=shell.summary.outlet[1].T - tubes.summary.outlet.T,
      absLevel=shell.phaseBorder.level_abs,
      relLevel=shell.phaseBorder.level_rel)) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-50,-92})));
protected
   ClaRa.Basics.Interfaces.EyeIn eye_int2 annotation (Placement(
        transformation(extent={{-51,-43},{-49,-41}})));
public
   ClaRa.Basics.Interfaces.EyeOut eye2 if showData annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-100,-40}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-110,0})));
protected
   ClaRa.Basics.Interfaces.EyeIn eye_int1
    annotation (Placement(transformation(extent={{27,-59},{29,-57}})));
public
   ClaRa.Basics.Interfaces.EyeOut eye1 if showData annotation (
      Placement(transformation(
        extent={{100,10},{120,30}},
        rotation=270,
        origin={20,0}), iconTransformation(
        extent={{100,10},{120,30}},
        rotation=270,
        origin={20,0})));
  ClaRa.Basics.Interfaces.FluidPortIn aux1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-110,70},{-90,90}})));
  ClaRa.Basics.Interfaces.FluidPortIn aux2(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-110,50},{-90,70}}),
        iconTransformation(extent={{-110,50},{-90,70}})));
  Modelica.Blocks.Interfaces.RealOutput level(value = if outputAbs then shell.summary.outline.level_abs else shell.summary.outline.level_rel) if levelOutput annotation (Placement(transformation(extent={{80,-90},{100,-110}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={80,-110})));
equation
   assert(diameter_o > diameter_i,
     "Outer diameter of tubes must be greater than inner diameter");

  connect(tubes.inlet, In2) annotation (Line(
      points={{62,20},{100,20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(tubes.outlet, Out2) annotation (Line(
      points={{42,20},{-100,20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));

   eye_int1.m_flow=-shell.outlet[1].m_flow;
   eye_int1.T=shell.summary.outlet[1].T-273.15;
   eye_int1.s=shell.fluidOut[1].s/1000;
   eye_int1.h=shell.summary.outlet[1].h/1000;
   eye_int1.p=shell.summary.outlet[1].p/100000;

   eye_int2.m_flow=-tubes.outlet.m_flow;
   eye_int2.T=tubes.summary.outlet.T-273.15;
   eye_int2.s=tubes.fluidOut.s/1000;
   eye_int2.h=tubes.summary.outlet.h/1000;
   eye_int2.p=tubes.summary.outlet.p/100000;

  connect(In1, shell.inlet[1]) annotation (Line(
      points={{0,98},{0,56},{1.77636e-015,56}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shell.outlet[1], Out1) annotation (Line(
      points={{-1.77636e-015,36},{-1.77636e-015,-26},{0,-26},{0,-100}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(aux1, shell.inlet[2]) annotation (Line(
      points={{-100,80},{1.77636e-015,80},{1.77636e-015,56}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(aux2, shell.inlet[3]) annotation (Line(
      points={{-100,60},{1.77636e-015,60},{1.77636e-015,56}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));

  connect(shell.heat[1], wall.outerPhase) annotation (Line(
      points={{9.5,46},{22.8667,46},{22.8667,45}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shell.heat[2], wall.outerPhase) annotation (Line(
      points={{10.5,46},{14,46},{14,45},{22.8667,45}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(tubes.heat, wall.innerPhase) annotation (Line(
      points={{52,30},{52,30},{52,44.8},{42.6,44.8}},
      color={167,25,48},
      thickness=0.5));
  connect(eye_int2, eye2) annotation (Line(points={{-50,-42},{-100,-42},{-100,-40}}, color={190,190,190}));
  connect(eye_int1, eye1) annotation (Line(points={{28,-58},{28,-58},{28,-110},{40,-110}},
                                                                                         color={190,190,190}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
                              Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}})));
end HEXvle2vle_L3_2ph_CH_simple;
