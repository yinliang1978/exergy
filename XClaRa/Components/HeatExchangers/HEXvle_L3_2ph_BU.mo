within Exergy.XClaRa.Components.HeatExchangers;
model HEXvle_L3_2ph_BU
  "Single side: VLE | L3 | two phase at shell side | Block shape | U-type"
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
  extends ClaRa.Basics.Icons.HEX04;
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=if not heatFlowIsLoss then -heat.Q_flow else 0,
    powerAux=0) if                                                                                                     contributeToCycleSummary;

  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");

  outer ClaRa.SimCenter simCenter;

  model Outline
    extends ClaRa.Basics.Icons.RecordIcon;
    parameter Boolean showExpertSummary=false;
    input ClaRa.Basics.Units.HeatFlowRate Q_flow "Heat flow rate";
    input ClaRa.Basics.Units.Length absLevel "Absolute filling level";
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

  //________________________________ Shell geometry _______________________________//
  parameter ClaRa.Basics.Units.Length
                      length=10 "Length of the HEX" annotation (Dialog(
      tab="Shell Side",
      group="Geometry", groupImage="modelica://ClaRa/figures/ParameterDialog/HEX_ParameterDialog_BUshell2ph2.png"));
  parameter ClaRa.Basics.Units.Length
                      height=3 "Height of HEX"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length
                      width=3 "Width of HEX"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length
                      z_in_shell=length/2 "Inlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
    parameter ClaRa.Basics.Units.Length
                        z_in_aux1=length/2
    "Inlet position of auxilliary1 from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length
                      z_in_aux2=length/2
    "Inlet position of auxilliary2 from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length
                      z_out_shell=length/2 "Outlet position from bottom"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length radius_flange=0.05
    "Flange radius of all flanges"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Mass
                    mass_struc=0
    "Mass of inner structure elements, additional to the tubes itself"
    annotation (Dialog(tab="Shell Side", group="Geometry"));

  //________________________________ Shell nominal parameter _____________________________________//
  parameter ClaRa.Basics.Units.MassFlowRate
                            m_flow_nom_shell=10
    "Nominal mass flow on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));
  parameter ClaRa.Basics.Units.Pressure
                        p_nom_shell=10 "Nominal pressure on shell side"
    annotation (Dialog(tab="Shell Side", group="Nominal Values"));

  //________________________________ Shell initialisation  _______________________________________//
   parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_liq_start=-10 +
      TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(
      medium_shell, p_start_shell) "Start specific enthalpy of liquid phase"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));
   parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_vap_start=+10 +
      TILMedia.VLEFluidFunctions.dewSpecificEnthalpy_pxi(medium_shell,
      p_start_shell) "Start specific enthalpy of steam phase"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));

  parameter ClaRa.Basics.Units.Pressure
                        p_start_shell=1e5 "Start value of sytsem pressure"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter Real level_rel_start=0.5 "Start value for relative filling Level" annotation (Dialog(tab="Shell Side", group="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initTypeShell=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation"
    annotation (Dialog(tab="Shell Side", group="Initialisation"));
  //________________________________ Shell epert settings  _______________________________________//
  parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
    "Smoothness of level calculation (table based)"                                                                                                     annotation (Dialog(tab="Shell Side", group="Expert Settings"));

  //*********************************** / TUBE SIDE \ ***********************************//
  //________________________________ Tubes geometry _______________________________//
//   parameter SI.Length d_i=0.048 "Inner diameter of horizontal tubes"
//     annotation (Dialog(
//       tab="Tubes",
//       group="Geometry",
//       groupImage="modelica://ClaRa/figures/ParameterDialog/HollowBlockWithTubes_2.png"));
  parameter ClaRa.Basics.Units.Length
                      diameter_o=0.05 "Outer diameter of horizontal tubes"
    annotation (Dialog(tab="Tubes", group="Geometry",
      groupImage="modelica://ClaRa/figures/ParameterDialog/HollowBlockWithTubes_2.png"));
  parameter Integer N_tubes=1000 "Number of horizontal tubes"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter Integer N_passes=1 "Number of passes of the internal tubes"
    annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter Boolean parallelTubes=false
    "True, if tubes are parallel to main orientation, else false"
    annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter Boolean staggeredAlignment=true
    "True, if the tubes are aligned staggeredly, false otherwise"
    annotation (Dialog(tab="Tubes", group="Geometry"));

  parameter ClaRa.Basics.Units.Length
                      Delta_z_par=2*diameter_o
    "Distance between tubes parallel to flow direction (center to center)"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter ClaRa.Basics.Units.Length
                      Delta_z_ort=2*diameter_o
    "Distance between tubes orthogonal to flow direction (center to center)"
    annotation (Dialog(tab="Tubes", group="Geometry"));
  parameter Integer N_rows=integer(ceil(sqrt(N_tubes))*N_passes)
    "Number of pipe rows in flow direction"                                                              annotation(Dialog(tab="Tubes", group="Geometry"));

  parameter Real CF_geo=1 "Correction coefficient due to fins etc."
    annotation (Dialog(tab="Tubes", group="Geometry"));

  //*********************************** / HOTWELL \ ***********************************//
  parameter ClaRa.Basics.Units.Length height_hotwell=1 "Height of the hotwell"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length width_hotwell=1 "Width of the hotwell"
    annotation (Dialog(tab="Shell Side", group="Geometry"));
  parameter ClaRa.Basics.Units.Length length_hotwell=1 "Length of the hotwell"
    annotation (Dialog(tab="Shell Side", group="Geometry"));

//*********************************** / EXPERT Settings and Visualisation \ ***********************************//
  parameter ClaRa.Basics.Units.Time Tau_cond=0.3
    "Time constant of condensation" annotation (Dialog(tab=
          "Expert Settings", group="Zone Interaction at Shell Side"));
  parameter ClaRa.Basics.Units.Time Tau_evap=0.03
    "Time constant of evaporation" annotation (Dialog(tab=
          "Expert Settings", group="Zone Interaction at Shell Side"));
  parameter ClaRa.Basics.Units.CoefficientOfHeatTransfer alpha_ph=50000
    "HTC of the phase border" annotation (Dialog(tab="Expert Settings",
        group="Zone Interaction at Shell Side"));
  parameter ClaRa.Basics.Units.Area A_phaseBorder=shell.geo.A_hor*100
    "Heat transfer area at phase border" annotation (Dialog(tab=
          "Expert Settings", group="Zone Interaction at Shell Side"));
  parameter Real expHT_phases=0
    "Exponent for volume dependency on inter phase HT"                             annotation (Dialog(tab="Expert Settings", group="Zone Interaction at Shell Side"));
  parameter Real absorbInflow=1
    "Absorption of incoming mass flow to the zones 1: perfect in the allocated zone, 0: perfect according to steam quality"
                                                                                              annotation (Dialog(tab="Expert Settings", group="Zone Interaction at Shell Side"));
  parameter Boolean equalPressures=true
    "True if pressure in liquid and vapour phase is equal"                                     annotation (Dialog(tab="Expert Settings", group="Zone Interaction at Shell Side"));

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "True, if expert summary should be applied"                                                               annotation (Dialog(tab="Summary and Visualisation"));
  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                              annotation (Dialog(tab="Summary and Visualisation"));
  parameter Boolean levelOutput = false
    "True, if Real level connector shall be addded"                                      annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean outputAbs = false "True, if absolute level is at output"  annotation(Dialog(enable = levelOutput, tab="Summary and Visualisation"));

  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean heatFlowIsLoss = true
    "True if heat flow is a loss (not a process product)"                                       annotation(Dialog(tab="Summary and Visualisation"));

  ClaRa.Basics.Interfaces.FluidPortOut outlet(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}})));
  ClaRa.Basics.Interfaces.FluidPortIn inlet(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-10,90},{10,110}}),
        iconTransformation(extent={{-10,90},{10,110}})));

  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_L3_TwoZonesNPort shell(
    medium=medium_shell,
    p_nom=p_nom_shell,
    redeclare model HeatTransfer = HeatTransfer_Shell,
    redeclare model PressureLoss = PressureLossShell,
    m_flow_nom=m_flow_nom_shell,
    useHomotopy=useHomotopy,
    p_start=p_start_shell,
    initType=initTypeShell,
    showExpertSummary=showExpertSummary,
    Tau_cond=Tau_cond,
    Tau_evap=Tau_evap,
    alpha_ph=alpha_ph,
    A_heat_ph=A_phaseBorder,
    h_liq_start=h_liq_start,
    h_vap_start=h_vap_start,
    level_rel_start=level_rel_start,
    exp_HT_phases=expHT_phases,
    heatSurfaceAlloc=2,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.RealSeparated
        (
        level_rel_start=level_rel_start,
        radius_flange=radius_flange,
        absorbInflow=absorbInflow,
        smoothness=smoothness),
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowBlockWithTubesAndHotwell
        (
        height=height,
        width=width,
        length=length,
        diameter_t=diameter_o,
        staggeredAlignment=staggeredAlignment,
        height_hotwell=height_hotwell,
        width_hotwell=width_hotwell,
        length_hotwell=length_hotwell,
        Delta_z_par=Delta_z_par,
        Delta_z_ort=Delta_z_ort,
        N_inlet=3,
        z_out={z_out_shell},
        flowOrientation=ClaRa.Basics.Choices.GeometryOrientation.vertical,
        z_in={z_in_shell,z_in_aux1,z_in_aux2},
        N_tubes=N_tubes,
        N_passes=N_passes,
        parallelTubes=parallelTubes,
        CF_geo={1,CF_geo},
        N_rows=N_rows),
    equalPressures=equalPressures) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={0,0})));

protected
TILMedia.VLEFluid_ph fluidIn(
    each vleFluidType=medium_shell,
    h=actualStream(inlet.h_outflow),
    p=inlet.p)                                                           annotation (Placement(transformation(extent={{-90,-12},
            {-70,8}},                                                                                                   rotation=0)));
TILMedia.VLEFluid_ph fluidOut(
    each vleFluidType=medium_shell,
    h=actualStream(outlet.h_outflow),
    p=outlet.p)                                                         annotation (Placement(transformation(extent={{-10,-70},
            {10,-50}},                                                                                                  rotation=0)));

public
    outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";

  Exergy.Utilities.ViewObjectNE viewObject(nEnergy={2,1,0,1});

  Exergy.Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));
equation

  //  viewObject.h[1].E_flow = inlet.m_flow*noEvent(actualStream(inlet.h_outflow));
  viewObject.h[1].E_flow = inlet.m_flow*fluidIn.h;
  viewObject.h[1].Ex_flow = inlet.m_flow*(fluidIn.h-refEnv.T*fluidIn.s);

  viewObject.h[2].E_flow = outlet.m_flow*fluidOut.h;
  viewObject.h[2].Ex_flow = outlet.m_flow*(fluidOut.h-refEnv.T*fluidOut.s);

    viewObject.q[1].E_flow = heat.Q_flow;
  viewObject.q[1].Ex_flow = heat.Q_flow*(1-refEnv.T/heat.T);

    viewObject.E[1].E = shell.summary.fluid.mass * shell.summary.fluid.h;
   viewObject.E[1].Ex = shell.summary.fluid.mass * (shell.summary.fluid.h - refEnv.T*{shell.liq.s,shell.vap.s});
   //viewObject.E[1].Ex = shell.summary.fluid.mass * ( refEnv.T*(shell.summary.fluid.rho));

  connect(viewObject.viewOutput,viewOutput);

public
  Summary summary(outline(
      showExpertSummary=showExpertSummary,
      Q_flow=sum(shell.heat.Q_flow),
      absLevel=shell.phaseBorder.level_abs,
      relLevel=shell.phaseBorder.level_rel)) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-50,-92})));
public
   ClaRa.Basics.Interfaces.EyeOut eye if showData annotation (Placement(
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
    annotation (Placement(transformation(extent={{-110,70},{-90,90}}),
        iconTransformation(extent={{-110,70},{-90,90}})));
  ClaRa.Basics.Interfaces.FluidPortIn aux2(Medium=medium_shell)
    annotation (Placement(transformation(extent={{-110,50},{-90,70}}),
        iconTransformation(extent={{-110,50},{-90,70}})));
  ClaRa.Basics.Interfaces.HeatPort_a heat
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  Adapters.Scalar2VectorHeatPort scalar2VectorHeatPort(N=2, final equalityMode="Equal Temperatures")
    annotation (Placement(transformation(extent={{60,-10},{40,10}})));

  Modelica.Blocks.Interfaces.RealOutput level(value = if outputAbs then shell.summary.outline.level_abs else shell.summary.outline.level_rel) if levelOutput annotation (Placement(transformation(extent={{78,-90},{98,-110}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={80,-110})));
equation
//   assert(diameter_o > d_i,
//     "Outer diameter of tubes must be greater than inner diameter");

   eye_int1.m_flow = shell.summary.outlet[1].m_flow;
   eye_int1.T = shell.summary.outlet[1].T - 273.15;
   eye_int1.s = shell.fluidOut[1].s/1e3;
   eye_int1.p = shell.outlet[1].p/1e5;
   eye_int1.h = shell.summary.outlet[1].h/1e3;

  connect(inlet, shell.inlet[1])
                               annotation (Line(
      points={{0,100},{0,10},{1.77636e-015,10}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(shell.outlet[1], outlet)
                                 annotation (Line(
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
  connect(scalar2VectorHeatPort.heatVector, shell.heat) annotation (Line(
      points={{40,0},{18,0},{18,-1.77636e-015},{10,-1.77636e-015}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(scalar2VectorHeatPort.heatScalar, heat) annotation (Line(
      points={{60,0},{100,0}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_int1, eye) annotation (Line(points={{28,-58},{28,-58},{28,-98},{28,-98}}, color={190,190,190}));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})),
                              Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}})));
end HEXvle_L3_2ph_BU;
