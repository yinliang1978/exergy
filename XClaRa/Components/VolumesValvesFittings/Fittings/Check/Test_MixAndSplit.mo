within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Check;
model Test_MixAndSplit
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

extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb60;

  Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Join_L2_flex
    flxibleJoin(
    h_nom=1900e3,
    p_nom(displayUnit="Pa") = 5000,
    h_start=1800e3,
    initType=ClaRa.Basics.Choices.Init.steadyState,
    N_ports_in=2,
    m_flow_in_nom={25,100},
    volume=0.05,
    preciseTwoPhase=false,
    p_start=3000000,
    showExpertSummary=true) annotation (Placement(transformation(
          extent={{-60,20},{-80,40}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_phxi massFlowSource_XRG(h_const=
        800e3, p_const=3000000)
    annotation (Placement(transformation(extent={{20,20},{0,40}})));
  inner ClaRa.SimCenter simCenter(redeclare replaceable
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1)
    annotation (Placement(transformation(extent={{80,158},{100,178}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource_XRG2(
    m_flow_const=43.551,
    variable_m_flow=true,
    h_const=1800e3) annotation (Placement(transformation(extent={{-120,
            20},{-100,40}})));
  Modelica.Blocks.Sources.Ramp ramp1(
    duration=200,
    height=200,
    offset=-150,
    startTime=200)
    annotation (Placement(transformation(extent={{-192,-94},{-172,-74}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource_XRG1(
    m_flow_const=43.551,
    variable_m_flow=true,
    h_const=3000e3,
    variable_h=true)
    annotation (Placement(transformation(extent={{20,50},{0,70}})));
  Modelica.Blocks.Sources.Ramp ramp2(
    duration=100,
    height=200,
    startTime=600,
    offset=-100)
    annotation (Placement(transformation(extent={{88,108},{68,128}})));
  ClaRa.Visualisation.Quadruple quadruple annotation (Placement(
        transformation(extent={{-84,2},{-124,12}})));
  Modelica.Blocks.Sources.Ramp ramp3(
    duration=100,
    startTime=1000,
    height=-2500e3,
    offset=3000e3)
    annotation (Placement(transformation(extent={{90,80},{70,100}})));
  Exergy.XClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1
    valveLinear_1_XRG(redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (
        m_flow_nom=(50),
        Delta_p_nom=(100000),
        rho_in_nom=if ((0) <> 0) then (0) else 10,
        Delta_p_eps=if ((0) > 0) then (0) else 100)) annotation (
      Placement(transformation(extent={{-40,24},{-20,36}})));
  Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Join_L2_Y
    join_Y(
    m_flow_in_nom={25,100},
    p_nom(displayUnit="Pa") = 5000,
    h_nom=1900e3,
    h_start=1800e3,
    initType=ClaRa.Basics.Choices.Init.steadyState,
    volume=0.05,
    preciseTwoPhase=false,
    p_start=3000000,
    showExpertSummary=true) annotation (Placement(transformation(
        extent={{-10,10},{10,-10}},
        rotation=180,
        origin={-70,-30})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource_XRG3(
    m_flow_const=43.551,
    variable_m_flow=true,
    h_const=1800e3) annotation (Placement(transformation(extent={{-120,
            -40},{-100,-20}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource_XRG4(
    m_flow_const=43.551,
    variable_m_flow=true,
    h_const=3000e3,
    variable_h=true) annotation (Placement(transformation(extent={{
            20,-10},{0,10}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_phxi massFlowSource_XRG5(h_const=
        800e3, p_const=3000000) annotation (Placement(
        transformation(extent={{20,-40},{0,-20}})));
  Exergy.XClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1
    valveLinear_1_XRG1(redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (
        m_flow_nom=(50),
        Delta_p_nom=(100000),
        rho_in_nom=if ((0) <> 0) then (0) else 10,
        Delta_p_eps=if ((0) > 0) then (0) else 100)) annotation (
      Placement(transformation(extent={{-40,-36},{-20,-24}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource_XRG6(
    m_flow_const=43.551,
    variable_m_flow=true,
    h_const=1800e3) annotation (Placement(transformation(extent={{-120,
            -100},{-100,-80}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource_XRG7(
    m_flow_const=43.551,
    variable_m_flow=true,
    h_const=3000e3,
    variable_h=true) annotation (Placement(transformation(extent={{
            20,-70},{0,-50}})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_phxi massFlowSource_XRG8(h_const=
        800e3, p_const=3000000) annotation (Placement(
        transformation(extent={{20,-100},{0,-80}})));
  Exergy.XClaRa.Components.VolumesValvesFittings.Valves.ValveVLE_L1
    valveLinear_1_XRG2(redeclare model PressureLoss =
        ClaRa.Components.VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint
        (
        m_flow_nom=(50),
        Delta_p_nom=(100000),
        rho_in_nom=if ((0) <> 0) then (0) else 10,
        Delta_p_eps=if ((0) > 0) then (0) else 100)) annotation (
      Placement(transformation(extent={{-40,-96},{-20,-84}})));
  Exergy.XClaRa.Components.VolumesValvesFittings.Fittings.Split_L2_Y
    split_Y(
    p_nom(displayUnit="Pa") = 5000,
    h_nom=1900e3,
    h_start=1800e3,
    initType=ClaRa.Basics.Choices.Init.steadyState,
    volume=0.05,
    m_flow_out_nom={25,100},
    preciseTwoPhase=false,
    p_start=3000000,
    showExpertSummary=true) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-70,-90})));
  ClaRa.Visualisation.Quadruple quadruple1 annotation (Placement(
        transformation(extent={{-84,-72},{-124,-62}})));
  ClaRa.Visualisation.Quadruple quadruple2 annotation (Placement(
        transformation(extent={{-84,-60},{-124,-50}})));
equation
  connect(ramp1.y, massFlowSource_XRG2.m_flow) annotation (Line(
      points={{-171,-84},{-144,-84},{-144,36},{-122,36}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(ramp2.y, massFlowSource_XRG1.m_flow) annotation (Line(
      points={{67,118},{38,118},{38,66},{22,66}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp3.y, massFlowSource_XRG1.h) annotation (Line(
      points={{69,90},{52,90},{52,60},{22,60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(valveLinear_1_XRG.outlet, massFlowSource_XRG.steam_a) annotation (
      Line(
      points={{-20,30},{0,30}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flxibleJoin.inlet[1], valveLinear_1_XRG.inlet)
                                                  annotation (Line(
      points={{-60,29.5},{-39.95,29.5},{-39.95,30},{-40,30}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flxibleJoin.inlet[2], massFlowSource_XRG1.steam_a)
                                                      annotation (Line(
      points={{-60,30.5},{-42,30.5},{-42,60},{0,60}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flxibleJoin.outlet, massFlowSource_XRG2.steam_a)
                                                    annotation (Line(
      points={{-80,30},{-100,30}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp1.y, massFlowSource_XRG3.m_flow) annotation (Line(
      points={{-171,-84},{-144,-84},{-144,-24},{-122,-24}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(massFlowSource_XRG3.steam_a, join_Y.outlet) annotation (Line(
      points={{-100,-30},{-80,-30}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp2.y, massFlowSource_XRG4.m_flow) annotation (Line(
      points={{67,118},{38,118},{38,6},{22,6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp3.y,massFlowSource_XRG4. h) annotation (Line(
      points={{69,90},{52,90},{52,0},{22,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource_XRG4.steam_a, join_Y.inlet2) annotation (Line(
      points={{0,0},{-70,0},{-70,-20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(valveLinear_1_XRG1.outlet, massFlowSource_XRG5.steam_a) annotation (
      Line(
      points={{-20,-30},{0,-30}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(ramp1.y, massFlowSource_XRG6.m_flow) annotation (Line(
      points={{-171,-84},{-122,-84}},
      color={0,0,127},
      smooth=Smooth.Bezier));
  connect(ramp2.y, massFlowSource_XRG7.m_flow) annotation (Line(
      points={{67,118},{38,118},{38,-54},{22,-54}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp3.y,massFlowSource_XRG7. h) annotation (Line(
      points={{69,90},{52,90},{52,-60},{22,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(valveLinear_1_XRG2.outlet, massFlowSource_XRG8.steam_a) annotation (
      Line(
      points={{-20,-90},{0,-90}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(split_Y.inlet, valveLinear_1_XRG2.inlet) annotation (Line(
      points={{-60,-90},{-40,-90}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(massFlowSource_XRG7.steam_a,split_Y. outlet2) annotation (Line(
      points={{0,-60},{-70,-60},{-70,-80}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(split_Y.outlet1, massFlowSource_XRG6.steam_a) annotation (Line(
      points={{-80,-90},{-100,-90}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(join_Y.inlet1, valveLinear_1_XRG1.inlet) annotation (Line(
      points={{-60,-30},{-49,-30},{-49,-30},{-40,-30}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(split_Y.eye[1], quadruple1.eye) annotation (Line(
      points={{-80,-83.5},{-80,-67},{-84,-67}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(quadruple2.eye, split_Y.eye[2]) annotation (Line(
      points={{-84,-55},{-80,-55},{-80,-84.5}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flxibleJoin.eye, quadruple.eye) annotation (Line(
      points={{-80,22},{-82,22},{-82,7},{-84,7}},
      color={190,190,190},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-200,-100},
            {100,180}}), graphics={                Text(
          extent={{-198,180},{-34,168}},
          lineColor={0,128,0},
          textString="Tested 01. Mar. 2013 //FG ",
          horizontalAlignment=TextAlignment.Left),
                                Text(
          extent={{-198,174},{22,52}},
          lineColor={0,128,0},
          textString="_______________________________________________________________________
PURPOSE:
Show different join or mixing models and illustrate the differences. The model compares a 
flexible join supporting an arbitrary number of inlet mass flows and a Y-type join having
two inlets. Additionally, a Y-type split is applied at the bottom of the diagram. As the
models are capable for reverse flow (as shown by changing the mass flow diraction after 
300 s) the split elements comes to the same results as the join models. However, it is
recommended to consider the design flow direction in the model choice as it eases the
initialisation.
_______________________________________________________________________
LOOK AT:
compare the corresponding outlet connector variables in the summaries. An important 
effect is the phase change at the outlet port after  flow reversal. This is due to 
the fact that the two mass flow sources have different specific enthalpies referring
to vapour and liquid phase, respectively. To calculate this high-gradient transients
rapidly the Boolean parameter preciseTwoPhase in the expert settings dialog is set to 
false. This refers to a cut link between the dynamic energy and mass balance equations,
i.e. the term h*V*der(rho) is ommited. This is okay in most cases since the dynamics of
a power plant do not depend strongly on relative small join elements.
Please note, setting preciseTwoPhase to true may induce strong mass flow oscillations
which can be suppressed by adding additional pressure losses.

_______________________________________________________________________
",        fontSize=10,
          horizontalAlignment=TextAlignment.Left)}),
    experiment(StopTime=3000),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=true)));
end Test_MixAndSplit;
