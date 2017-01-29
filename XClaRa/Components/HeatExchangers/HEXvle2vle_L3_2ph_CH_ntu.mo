within Exergy.XClaRa.Components.HeatExchangers;
model HEXvle2vle_L3_2ph_CH_ntu
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
    input Real effectiveness[3] if showExpertSummary "Effectivenes of HEX";
    input Real kA[3](unit="W/K") if showExpertSummary "Overall heat resistance";
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
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching);
//, groupImage="modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialog_CHgeneral.png"
  //________________________________ Shell geometry _______________________________//
  parameter SI.Length length=10 "Length of the HEX" annotation (Dialog(
      tab="Shell Side",
      group="Geometry",
      groupImage=
          "modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialog_CHshell.png"));
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
  parameter SI.Mass mass_struc=0
    "Mass of inner structure elements, additional to the tubes itself"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Choices.GeometryOrientation orientation=ClaRa.Basics.Choices.GeometryOrientation.vertical
    "Orientation of the component"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Choices.GeometryOrientation flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.vertical
    "Flow orientation at shell side"
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
      p_start_shell) "Start value of liquid specific enthalpy" annotation (Dialog(tab="Shell Side", group="Initialisation"));
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
    "Medium to be used for tubes flow" annotation (choices(
      choice=simCenter.fluid1 "First fluid defined in global simCenter",
      choice=simCenter.fluid2 "Second fluid defined in global simCenter",
      choice=simCenter.fluid3 "Third fluid defined in global simCenter"),
                                                          Dialog(tab="Tubes",
        group="Fundamental Definitions"));
  replaceable model HeatTransferTubes =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.CharLine_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.TubeType_L2
    "Heat transfer mode at the tubes side" annotation (Dialog(tab="Tubes",
        group="Fundamental Definitions"), choicesAllMatching);
  replaceable model PressureLossTubes =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.TubeType_L2
    "Pressure loss model at the tubes side" annotation (Dialog(tab="Tubes",
        group="Fundamental Definitions"), choicesAllMatching);

  //________________________________ Tubes geometry _______________________________//

  parameter SI.Length diameter_i=0.048 "Inner diameter of horizontal tubes"
                                         annotation (Dialog(
      tab="Tubes",
      group="Geometry",
      groupImage=
          "modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialogTubes.png"));
  parameter SI.Length diameter_o=0.05 "Outer diameter of horizontal tubes"
    annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter SI.Length length_tubes=10 "Length of the tubes (one pass)"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter Integer N_tubes=1000 "Number of horizontal tubes" annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter Integer N_passes=1 "Number of passes of the internal tubes"
                                                                       annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter SI.Length z_in_tubes=length/2 "Inlet position from bottom"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter SI.Length z_out_tubes=length/2 "Outlet position from bottom"
    annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter Boolean parallelTubes=true
    "True, if tubes are parallel to shell flow orientation"                                    annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter Boolean staggeredAlignment=true
    "True, if the tubes are aligned staggeredly"                                         annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter SI.Length Delta_z_par=2*diameter_o
    "Distance between tubes parallel to flow direction (center to center)"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter SI.Length Delta_z_ort=2*diameter_o
    "Distance between tubes orthogonal to flow direction (center to center)"
    annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter Integer N_rows=integer(ceil(sqrt(N_tubes))*N_passes)
    "Number of pipe rows in shell flow direction"                                                              annotation(Dialog(tab="Tubes", group="Geometry"));
  parameter Real CF_geo=1 "Correction coefficient due to fins etc." annotation (Dialog(tab="Tubes", group="Geometry"));

  //________________________________ Tubes nominal parameter _____________________________________//
  parameter SI.MassFlowRate m_flow_nom_tubes=10
    "Nominal mass flow on tubes side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
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

  //________________________________ Wall initialisation _______________________________________//
  parameter SI.Temperature T_w_tube_start[3]=ones(3)*293.15
    "Initial temperature at inner phase"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));
  parameter SI.Temperature T_w_shell_start[3]=ones(3)*293.15
    "Initial temperature at outer phase"
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
    "Smoothness of level calculation (table based)"                                                                                                     annotation (Dialog(tab="Expert Settings", group="Mass Accumulation at Shell Side"));

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "True, if expert summary should be applied"                                                               annotation (Dialog(tab="Summary and Visualisation"));
  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                              annotation (Dialog(tab="Summary and Visualisation"));
  parameter Boolean levelOutput = false
    "True, if Real level connector shall be addded"                                      annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean outputAbs = false "True, if absolute level is at output"  annotation(Dialog(enable = levelOutput, tab="Summary and Visualisation"));

  parameter Real gain_eff=1
    "Avoid effectiveness > 1, high gain_eff leads to stricter observation but may cause numeric errors"
                                                                                              annotation (Dialog(tab="Expert Settings", group="NTU model"));
  parameter SI.Time Tau_stab=0.1
    "Time constant for numeric stabilisation w.r.t. heat flow rates"
    annotation (Dialog(tab="Expert Settings", group="NTU model"));
  replaceable function HeatCapacityAveraging =
      ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.Functions.ArithmeticMean
    constrainedby
    ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.Functions.GeneralMean
    "Method for Averaging of heat capacities"
    annotation (Dialog(tab="Expert Settings", group="NTU model"),choicesAllMatching);

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
        Basics.ControlVolumes.Fundamentals.Geometry.PipeGeometry (
        z_in={z_in_tubes},
        z_out={z_out_tubes},
        diameter=diameter_i,
        N_passes=N_passes,
        length=length_tubes,
        N_tubes=N_tubes))
    annotation (Placement(transformation(extent={{84,10},{64,30}})));

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
    A_heat_ph=A_phaseBorder,
    exp_HT_phases=expHT_phases,
    heatSurfaceAlloc=2,
    redeclare model PhaseBorder =
        Basics.ControlVolumes.Fundamentals.SpacialDistribution.RealSeparated (
        level_rel_start=level_rel_start,
        absorbInflow=absorbInflow,
        radius_flange=radius_flange,
        smoothness=smoothness),
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowCylinderWithTubes
        (
        z_out={z_out_shell},
        length=length,
        staggeredAlignment=staggeredAlignment,
        Delta_z_par=Delta_z_par,
        Delta_z_ort=Delta_z_ort,
        diameter=diameter,
        final N_passes=N_passes,
        length_tubes=length_tubes,
        diameter_t=diameter_o,
        N_inlet=3,
        z_in={z_in_shell,z_in_aux1,z_in_aux2},
        N_tubes=N_tubes,
        orientation=orientation,
        CF_geo={1,CF_geo},
        flowOrientation=flowOrientation,
        N_rows=N_rows,
        final parallelTubes=parallelTubes,
        final N_baffle=0),
    equalPressures=equalPressures)
                           annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,42})));

public
  Summary summary(outline(
      showExpertSummary=showExpertSummary,
      Q_flow=sum(shell.heat.Q_flow),
      Delta_T_in=shell.summary.inlet[1].T - tubes.summary.inlet.T,
      Delta_T_out=shell.summary.outlet[1].T - tubes.summary.outlet.T,
      effectiveness=wall.summary.effectiveness,
      kA=wall.summary.kA,
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
        origin={-100,-42}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-110,0})));
protected
   ClaRa.Basics.Interfaces.EyeIn eye_int1
    annotation (Placement(transformation(extent={{27,-59},{29,-57}})));
public
   ClaRa.Basics.Interfaces.EyeOut eye1 if showData annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={28,-98}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,-110})));
  ClaRa.Basics.Interfaces.FluidPortIn aux1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-110,70},{-90,90}})));
  ClaRa.Basics.Interfaces.FluidPortIn aux2(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-110,50},{-90,70}}),
        iconTransformation(extent={{-110,50},{-90,70}})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.NTU_L3_standalone wall(
    length=length,
    redeclare replaceable model Material = WallMaterial,
    N_p=N_passes,
    radius_i=diameter_i/2,
    m_flow_i=tubes.inlet.m_flow,
    alpha_i=ones(3)*tubes.heattransfer.alpha,
    mass_struc=mass_struc,
    radius_o=diameter_o/2,
    h_i_inlet=inStream(tubes.inlet.h_outflow),
    h_o_inlet=inStream(shell.inlet[1].h_outflow),
    m_flow_o=sum(shell.inlet.m_flow),
    redeclare function HeatCapacityAveraging = HeatCapacityAveraging,
    medium_shell=medium_shell,
    medium_tubes=medium_tubes,
    T_w_i_start=T_w_tube_start,
    T_w_o_start=T_w_shell_start,
    initChoice=initTypeWall,
    outerPhaseChange=true,
    p_i=tubes.outlet.p,
    showExpertSummary=showExpertSummary,
    gain_eff=gain_eff,
    N_t=N_tubes,
    redeclare replaceable model HeatExchangerType =
        ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.HeatExchangerTypes.CounterFlow_L3,
    alpha_o={shell.heattransfer.alpha[2],shell.heattransfer.alpha[2],
        shell.heattransfer.alpha[1]},
    Tau_stab=Tau_stab,
    CF_geo=CF_geo,
    p_o=shell.inlet[1].p) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={54,40})));

  Adapters.Scalar2VectorHeatPort reallocateHeatFlows(final equalityMode="Equal Temperatures")
    annotation (Placement(transformation(extent={{16,40},{36,60}})));
  Modelica.Blocks.Interfaces.RealOutput level(value = if outputAbs then shell.summary.outline.level_abs else shell.summary.outline.level_rel) if levelOutput annotation (Placement(transformation(extent={{204,-126},{224,-106}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={80,-110})));
equation
   assert(diameter_o > diameter_i,
     "Outer diameter of tubes must be greater than inner diameter");

initial equation
  //        wall.T=(tubes.bulk.T+shell.bulk.T)/2;

equation
  connect(tubes.inlet, In2) annotation (Line(
      points={{84,20},{100,20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(tubes.outlet, Out2) annotation (Line(
      points={{64,20},{-100,20}},
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
      points={{0,98},{0,52},{1.77636e-015,52}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shell.outlet[1], Out1) annotation (Line(
      points={{-1.77636e-015,32},{-1.77636e-015,-26},{0,-26},{0,-100}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(aux1, shell.inlet[2]) annotation (Line(
      points={{-100,80},{1.77636e-015,80},{1.77636e-015,52}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(aux2, shell.inlet[3]) annotation (Line(
      points={{-100,60},{1.77636e-015,60},{1.77636e-015,52}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));

  connect(wall.innerPhase[1], tubes.heat) annotation (Line(
      points={{63.6667,40},{63.6667,38},{74,38},{74,30}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(wall.innerPhase[2], tubes.heat) annotation (Line(
      points={{63,40},{63,34},{74,34},{74,30}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(wall.innerPhase[3], tubes.heat) annotation (Line(
      points={{62.3333,40},{74,40},{74,30}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_int2, eye2) annotation (Line(points={{-50,-42},{-100,-42},{-100,-42}}, color={190,190,190}));
  connect(eye_int1, eye1) annotation (Line(points={{28,-58},{28,-98},{28,-98}}, color={190,190,190}));
  connect(wall.outerPhase, reallocateHeatFlows.heatVector) annotation (Line(
      points={{45,40},{42,40},{42,50},{36,50}},
      color={167,25,48},
      thickness=0.5));
  connect(reallocateHeatFlows.heatScalar, shell.heat[1]) annotation (Line(
      points={{16,50},{12,50},{12,42},{9.5,42}},
      color={167,25,48},
      thickness=0.5));
  connect(reallocateHeatFlows.heatScalar, shell.heat[2]) annotation (Line(
      points={{16,50},{16,50},{16,42},{10.5,42}},
      color={167,25,48},
      thickness=0.5));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics={Text(
          extent={{-86,78},{86,38}},
          lineColor={27,36,42},
          textString="NTU")}),Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={Text(
          extent={{4,14},{54,-10}},
          lineColor={118,124,127},
          textString="Heat flow from subcooled zone calculated
 by NTU-model is mixed with the heat flow 
from the other zones and distributed according
to temperature differences to the shell volume.
This is done due to numerical reasons.
In HEXvle2vle_L3_2ph_BU_ntu this is 
handled differently.")}));
end HEXvle2vle_L3_2ph_CH_ntu;
