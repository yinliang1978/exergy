within Exergy.XClaRa.Components.HeatExchangers;
model HEXvle2vle_L3_2ph_CU_ntu
  "VLE 2 VLE | L3 | 2 phase at shell side | Cylinder shape | U-type | NTU ansatz"
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

  //pure cross current is assumed. however, the geometry is quite complex and the idea of moving boundary cell method with three zones at each side might be extended to get better results

  // Extends from... ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  extends ClaRa.Basics.Icons.HEX04;
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
  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching);

  //*********************************** / SHELL SIDE \ ***********************************//
  //________________________________ Shell fundamentals _______________________________//
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_shell=simCenter.fluid1
    "Medium to be used for shell side" annotation (choices(
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

  //________________________________ Shell geometry _______________________________//
  parameter SI.Length length=10 "Length of the HEX" annotation (Dialog(
      tab="Shell Side",
      group="Geometry",
      groupImage=
          "modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialog_CUshell.png"));
  parameter SI.Length diameter=3 "Diameter of HEX"
    annotation (Dialog(tab="Shell Side", group="Geometry"));

  parameter ClaRa.Basics.Choices.GeometryOrientation orientation=ClaRa.Basics.Choices.GeometryOrientation.vertical
    "Orientation of the component"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Choices.GeometryOrientation flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.vertical
    "Flow orientation at shell side"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  final parameter Boolean parallelTubes = if orientation == flowOrientation and N_baffle==0 then true elseif orientation <> flowOrientation and N_baffle==0 then true else false
    "True if tues are parallel to shell flow orientation";

  parameter SI.Length z_in_shell=length/2 "Inlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));

  parameter SI.Length z_in_aux1=length/2
    "Inlet position of auxilliary1 from bottom"                                      annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter SI.Length z_in_aux2=length/2
    "Inlet position of auxilliary2 from bottom"                                     annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter SI.Length z_out_shell=length/2 "Outlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter SI.Length radius_flange=0.05 "Flange radius" annotation (Dialog(tab="Shell Side", group="Geometry"));

  parameter Integer N_baffle=0 "Number of baffles at shell side" annotation (Dialog(tab="Shell Side", group="Geometry"));
  final parameter SI.Mass mass_struc=0
    "Mass of inner structure elements, additional to the tubes itself"
    annotation (Dialog(tab="Shell Side", group="Geometry"));

  //________________________________ Shell nominal parameter _____________________________________//
  parameter SI.MassFlowRate m_flow_nom_shell=10
    "Nominal mass flow at shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));
  parameter SI.Pressure p_nom_shell=10 "Nominal pressure at shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));
  parameter SI.EnthalpyMassSpecific h_nom_shell=100e3
    "Nominal specific enthalpy at shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));

  //________________________________ Shell initialisation  _______________________________________//

    parameter SI.EnthalpyMassSpecific h_liq_start=-10 +
      TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(medium_shell,
      p_start_shell) "Start value of liquid specific enthalpy"
                                                              annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter SI.EnthalpyMassSpecific h_vap_start=+10 +
      TILMedia.VLEFluidFunctions.dewSpecificEnthalpy_pxi(medium_shell, p_start_shell)
    "Start value of vapour specific enthalpy"                                                 annotation (Dialog(tab="Shell Side", group="Initialisation"));

  parameter SI.Pressure p_start_shell=1e5 "Start value of shell fluid pressure"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter Real level_rel_start=0.5 "Start value for relative filling Level" annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeShell=ClaRa.Basics.Choices.Init.noInit
    "Type of shell fluid initialisation"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));

  //*********************************** / TUBE SIDE \ ***********************************//
  //________________________________ Tubes fundamentals _______________________________//
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid medium_tubes=simCenter.fluid1
    "Medium to be used for tubes" annotation (choices(
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
  parameter ClaRa.Basics.Units.Length diameter_i=0.048
    "Inner diameter of horizontal tubes"                                                    annotation (Dialog(tab="Tubes", group="Geometry",groupImage="modelica://ClaRa/figures/ParameterDialog/HollowBlockWithTubes_2.png"));
  parameter ClaRa.Basics.Units.Length diameter_o=0.05
    "Outer diameter of horizontal tubes"                                                   annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length length_tubes=10
    "Length of the tubes (one pass)"                                                   annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter Integer N_tubes=1000 "Number of horizontal tubes" annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter Integer N_passes=1 "Number of passes of the internal tubes" annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length z_in_tubes=length/2
    "Inlet position from bottom"                                                       annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length z_out_tubes=length/2
    "Outlet position from bottom"                                                        annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter Boolean staggeredAlignment=true
    "True, if the tubes are aligned staggeredly, false otherwise"                                         annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length Delta_z_par=2*diameter_o
    "Distance between tubes parallel to flow direction (center to center)"                    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length Delta_z_ort=2*diameter_o
    "Distance between tubes orthogonal to flow direction (center to center)"                  annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter Integer N_rows=integer(ceil(sqrt(N_tubes))*N_passes)
    "Number of pipe rows in flow direction"                                                              annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter Real CF_geo=1 "Correction coefficient due to fins etc." annotation (Dialog(tab="Tubes", group="Geometry"));

  //________________________________ Tubes nominal parameter _____________________________________//
  parameter ClaRa.Basics.Units.MassFlowRate m_flow_nom_tubes=10
    "Nominal mass flow on tube side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter ClaRa.Basics.Units.Pressure p_nom_tubes=10
    "Nominal pressure on tube side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_nom_tubes=10
    "Nominal specific enthalpy on tube side"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));
  parameter ClaRa.Basics.Units.HeatFlowRate Q_flow_nom=1e6
    "Nominal heat flow rate"
    annotation (Dialog(tab="Tubes", group="Nominal Values"));

  //________________________________ Tubes initialisation _______________________________________//
  parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_start_tubes=1e5
    "Start value of tube fluid specific enthalpy"                                                                   annotation (Dialog(tab="Tubes", group="Initialisation"));
  parameter ClaRa.Basics.Units.Pressure p_start_tubes=1e5
    "Start value of tube fluid pressure"                                                       annotation (Dialog(tab="Tubes", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeTubes=ClaRa.Basics.Choices.Init.noInit
    "Type of tube fluid initialisation"                                                                                  annotation (Dialog(tab="Tubes", group="Initialisation"));

  //*********************************** / WALL \ ***********************************//
  //________________________________ Wall fundamentals _______________________________//
  replaceable model WallMaterial =
      TILMedia.SolidTypes.TILMedia_Aluminum
    constrainedby TILMedia.SolidTypes.BaseSolid "Material of the cylinder"
    annotation (choicesAllMatching=true, Dialog(tab="Tube Wall", group =  "Fundamental Definitions"));
  //________________________________ Wall initialisation _______________________________________//
  parameter SI.Temperature T_w_i_start[3]=ones(3)*293.15
    "Initial temperature at inner phase"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));
  parameter SI.Temperature T_w_o_start[3]=ones(3)*293.15
    "Initial temperature at outer phase"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeWall=ClaRa.Basics.Choices.Init.noInit
    "Init option of Tube wall"
    annotation (Dialog(tab="Tube Wall", group="Initialisation"));

  //*********************************** / EXPERT Settings and Visualisation \ ***********************************//
  parameter SI.Time Tau_cond=0.3
    "|Expert Settings|Inter-Phase Heat and Mass Transfer at Shell Side|Time constant of condensation";
  parameter SI.Time Tau_evap=0.03
    "|Expert Settings|Inter-Phase Heat and Mass Transfer at Shell Side|Time constant of evaporation";
  parameter SI.CoefficientOfHeatTransfer alpha_ph=50000
    "|Expert Settings|Inter-Phase Heat and Mass Transfer at Shell Side|HTC of the phase border";
  parameter SI.Area A_phaseBorder=shell.geo.A_hor*100
    "|Expert Settings|Inter-Phase Heat and Mass Transfer at Shell Side|Heat transfer area at phase border";
  parameter Real expHT_phases=0
    "|Expert Settings|Inter-Phase Heat and Mass Transfer at Shell Side|Exponent for volume dependency on inter phase HT";
  parameter Real absorbInflow=1
    "|Expert Settings|Inter-Phase Heat and Mass Transfer at Shell Side|Absorption of incoming mass flow to the zones 1: perfect in the allocated zone, 0: perfect according to steam quality";
  parameter Boolean equalPressures=true
    "True if pressure in liquid and vapour phase is equal"                                     annotation (Dialog(tab="Expert Settings", group="Zone Interaction at Shell Side"));

  parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
    "Smoothness of level calculation (table based)"                                                                                                     annotation (Dialog(tab="Expert Settings", group="Mass Accumulation at Shell Side"));

  replaceable function HeatCapacityAveraging =
      ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.Functions.ArithmeticMean
    constrainedby
    ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.Functions.GeneralMean
    "|Expert Settings|NTU model|Method for Averaging of heat capacities"
    annotation (choicesAllMatching);
  parameter Real gain_eff=1
    "|Expert Settings|NTU model|Avoid effectiveness > 1, high gain_eff leads to stricter observation but may cause numeric errors";
  parameter SI.Time Tau_stab=0.1
    "|Expert Settings|NTU model|Time constant for numeric stabilisation w.r.t. heat flow rates";

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "|Summary and Visualisation||True, if expert summary should be applied";
  parameter Boolean showData=true
    "|Summary and Visualisation||True, if a data port containing p,T,h,s,m_flow shall be shown, else false";
  parameter Boolean levelOutput = false
    "True, if Real level connector shall be addded"                                      annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean outputAbs = false "True, if absolute level is at output"  annotation(Dialog(enable = levelOutput, tab="Summary and Visualisation"));
  ClaRa.Basics.Interfaces.FluidPortIn In2(Medium=medium_tubes)
    annotation (Placement(transformation(extent={{90,-50},{110,-30}}),
        iconTransformation(extent={{90,-50},{110,-30}})));
  ClaRa.Basics.Interfaces.FluidPortOut Out2(Medium=medium_tubes)
    annotation (Placement(transformation(extent={{90,50},{110,70}}),
        iconTransformation(extent={{90,50},{110,70}})));
  ClaRa.Basics.Interfaces.FluidPortOut Out1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
  ClaRa.Basics.Interfaces.FluidPortIn In1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,88},{10,108}})));

  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_2 tubes(
    final medium=medium_tubes,
    final p_nom=p_nom_tubes,
    final h_nom=h_nom_tubes,
    final m_flow_nom=m_flow_nom_tubes,
    final useHomotopy=useHomotopy,
    final h_start=h_start_tubes,
    final p_start=p_start_tubes,
    final initType=initTypeTubes,
    redeclare model HeatTransfer = HeatTransferTubes,
    redeclare model PressureLoss = PressureLossTubes,
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.PipeGeometry (
        z_in={z_in_tubes},
        z_out={z_out_tubes},
        diameter=diameter_i,
        N_tubes=N_tubes,
        N_passes=N_passes,
        length=length_tubes),
    showExpertSummary=showExpertSummary,
    redeclare model PhaseBorder =
        Basics.ControlVolumes.Fundamentals.SpacialDistribution.IdeallyStirred)
    annotation (Placement(transformation(
        extent={{10,10},{-10,-10}},
        rotation=270,
        origin={84,0})));

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
        final radius_flange=radius_flange,
        final absorbInflow=absorbInflow,
        final smoothness=smoothness),
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowCylinderWithTubes
        (
        final z_out={z_out_shell},
        final length=length,
        final N_tubes=N_tubes,
        final z_in={z_in_shell,z_in_aux1,z_in_aux2},
        final staggeredAlignment=staggeredAlignment,
        final Delta_z_par=Delta_z_par,
        final Delta_z_ort=Delta_z_ort,
        final diameter=diameter,
        final N_passes=N_passes,
        final length_tubes=length_tubes,
        final diameter_t=diameter_o,
        final N_baffle=N_baffle,
        final N_rows=N_rows,
        final N_inlet=3,
        final orientation=orientation,
        final parallelTubes=parallelTubes,
        final CF_geo={1,CF_geo},
        final flowOrientation=flowOrientation),
    equalPressures=equalPressures) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,0})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.NTU_L3_standalone wall(
    redeclare replaceable model Material = WallMaterial,
    N_p=N_passes,
    radius_i=diameter_i/2,
    m_flow_i=tubes.inlet.m_flow,
    mass_struc=mass_struc,
    N_t=N_tubes,
    radius_o=diameter_o/2,
    redeclare function HeatCapacityAveraging = HeatCapacityAveraging,
    medium_shell=medium_shell,
    medium_tubes=medium_tubes,
    T_w_i_start=T_w_i_start,
    T_w_o_start=T_w_o_start,
    initChoice=initTypeWall,
    outerPhaseChange=true,
    alpha_o={shell.heattransfer.alpha[2],shell.heattransfer.alpha[2],shell.heattransfer.alpha[1]},
    alpha_i={1,1,1}*tubes.heattransfer.alpha,
    h_i_inlet=tubes.fluidIn.h,
    p_i=(tubes.fluidIn.p + tubes.fluidOut.p)/2,
    gain_eff=gain_eff,
    length=length_tubes,
    h_o_inlet=shell.fluidIn[1].h,
    m_flow_o=shell.inlet[1].m_flow,
    showExpertSummary=showExpertSummary,
    redeclare model HeatExchangerType =
        ClaRa.Basics.ControlVolumes.SolidVolumes.Fundamentals.HeatExchangerTypes.CrossFlow_L3,
    Tau_stab=Tau_stab,
    p_o=shell.inlet[1].p)
    "{shell.heattransfer.alpha[2],shell.heattransfer.alpha[2],shell.heattransfer.alpha[1]}"
                                                                                              annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=90,
        origin={57,0})));

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
  ClaRa.Basics.Interfaces.EyeOut eye2 if showData annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-100,-42}), iconTransformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={110,80})));
public
  ClaRa.Basics.Interfaces.EyeOut eye1 if showData annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={28,-98}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={40,-110})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int1
    annotation (Placement(transformation(extent={{27,-59},{29,-57}})));

public
  ClaRa.Basics.Interfaces.FluidPortIn aux1(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-110,70},{-90,90}})));
  ClaRa.Basics.Interfaces.FluidPortIn aux2(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-110,50},{-90,70}}),
        iconTransformation(extent={{-110,50},{-90,70}})));
  Adapters.Scalar2VectorHeatPort reallocateHeatFlows(final equalityMode="Equal Temperatures")
    annotation (Placement(transformation(extent={{20,-4},{40,16}})));

  Modelica.Blocks.Interfaces.RealOutput level(value = if outputAbs then shell.summary.outline.level_abs else shell.summary.outline.level_rel) if levelOutput annotation (Placement(transformation(extent={{204,-126},{224,-106}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={80,-110})));
equation
  assert(diameter_o > diameter_i,
    "Outer diameter of tubes must be greater than inner diameter");

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

  connect(tubes.inlet, In2) annotation (Line(
      points={{84,-10},{84,-40},{100,-40}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(tubes.outlet, Out2) annotation (Line(
      points={{84,10},{84,60},{100,60}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(wall.innerPhase[1], tubes.heat) annotation (Line(
      points={{66.6667,-8.88178e-016},{66.6667,-1.77636e-015},{74,-1.77636e-015}},
      color={191,0,0},
      smooth=Smooth.None,
      thickness=0.5));
  connect(wall.innerPhase[2], tubes.heat) annotation (Line(
      points={{66,-4.44089e-016},{66,-1.77636e-015},{74,-1.77636e-015}},
      color={191,0,0},
      smooth=Smooth.None,
      thickness=0.5));
  connect(wall.innerPhase[3], tubes.heat) annotation (Line(
      points={{65.3333,-4.44089e-016},{65.3333,1.77636e-015},{74,1.77636e-015}},
      color={191,0,0},
      smooth=Smooth.None,
      thickness=0.5));
  connect(eye_int1, eye1) annotation (Line(
      points={{28,-58},{28,-98}},
      color={190,190,190},
      smooth=Smooth.None));

  connect(In1, shell.inlet[1]) annotation (Line(
      points={{0,98},{0,10},{1.77636e-015,10}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shell.outlet[1], Out1) annotation (Line(
      points={{-1.77636e-015,-10},{-1.77636e-015,-30},{0,-30},{0,-100}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(aux1, shell.inlet[2]) annotation (Line(
      points={{-100,80},{1.77636e-015,80},{1.77636e-015,10}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(aux2, shell.inlet[3]) annotation (Line(
      points={{-100,60},{1.77636e-015,60},{1.77636e-015,10}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(reallocateHeatFlows.heatVector, wall.outerPhase) annotation (Line(
      points={{40,6},{44,6},{44,0},{48,0},{48,6.66134e-016}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(reallocateHeatFlows.heatScalar, shell.heat[1]) annotation (Line(
      points={{20,6},{14,6},{14,-1.77636e-015},{9.5,-1.77636e-015}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(reallocateHeatFlows.heatScalar, shell.heat[2]) annotation (Line(
      points={{20,6},{18,6},{18,-1.77636e-015},{10.5,-1.77636e-015}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_int2, eye2) annotation (Line(points={{-50,-42},{-100,-42},{-100,-42}}, color={190,190,190}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,
            -100},{100,100}}),
                         graphics={      Text(
          extent={{-24,64},{92,20}},
          lineColor={27,36,42},
          textString="NTU")}), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={Text(
          extent={{12,42},{62,18}},
          lineColor={118,124,127},
          textString="Heat flow from subcooled zone calculated
 by NTU-model is mixed with the heat flow 
from the other zones and distributed according
to temperature differences to the shell volume.
This is done due to numerical reasons.
In HEXvle2vle_L3_2ph_BU_ntu this is 
handled differently.")}));
end HEXvle2vle_L3_2ph_CU_ntu;
