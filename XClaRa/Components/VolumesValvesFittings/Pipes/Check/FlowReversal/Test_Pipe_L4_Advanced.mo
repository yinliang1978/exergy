within Exergy.XClaRa.Components.VolumesValvesFittings.Pipes.Check.FlowReversal;
model Test_Pipe_L4_Advanced
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

  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb50;

  Modelica.Blocks.Math.MultiSum multiSum(nu=2) annotation (Placement(
        transformation(
        extent={{-6,-6},{6,6}},
        rotation=180,
        origin={77,-24})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource(
    m_flow_const=0.1,
    variable_m_flow=true,
    h_const=200e3,
    m_flow_nom=0,
    variable_h=true,
    p_nom=1000) annotation (Placement(transformation(extent={{60,
            -39},{40,-19}})));
  inner ClaRa.SimCenter simCenter(redeclare replaceable
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1,
      useHomotopy=true) annotation (Placement(transformation(
          extent={{-80,-106},{-60,-86}})));
  PipeFlowVLE_L4_Advanced tube(
    z_in=0,
    z_out=0,
    showExpertSummary=true,
    showData=true,
    m_flow_nom=100,
    diameter_i=0.5,
    Delta_x=ones(tube.N_cv)*tube.length/tube.N_cv,
    h_start=ones(tube.N_cv)*200e3,
    N_cv=50,
    length=50,
    redeclare model PressureLoss =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L4,
    p_start=linspace(
        320000,
        300000,
        tube.N_cv),
    frictionAtInlet=true,
    suppressHighFrequencyOscillations=true,
    initType=ClaRa.Basics.Choices.Init.steadyState,
    frictionAtOutlet=true) annotation (Placement(transformation(extent={{20,-36},{-10,-24}})));

  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_phxi massFlowSink(
    variable_p=true,
    m_flow_nom=100,
    p_const=1000000,
    h_const=200e3,
    Delta_p=100000) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-46,-29})));
  inner Modelica.Fluid.System system
    annotation (Placement(transformation(extent={{60,-104},{80,-84}})));
  Modelica.Blocks.Sources.Step inlet_pressure(
    offset=1e5,
    startTime=100,
    height=1e4)
    annotation (Placement(transformation(extent={{-92,-45},{-72,-25}})));
  Modelica.Blocks.Sources.Ramp mass_flow_1(
    duration=1,
    height=10,
    offset=100,
    startTime=500) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={50,38})));

  Modelica.Blocks.Sources.Ramp T_wall(
    offset=293.15,
    startTime=1000,
    duration=100,
    height=30) annotation (Placement(transformation(
        extent={{-10.5,-10.5},{10.5,10.5}},
        rotation=0,
        origin={-81.5,7.5})));
  Utilities.Blocks.RealInputMultiplyer realInputMultiplyer(N=tube.N_cv) annotation (Placement(transformation(extent={{-46,-2},{-32,17}})));

  Modelica.Blocks.Sources.Ramp mass_flow_2(
    offset=100,
    startTime=1500,
    height=-300,
    duration=200) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={50,6})));

  Modelica.Blocks.Sources.Step inlet_pressure1(
    height=20e3,
    offset=200e3,
    startTime=200) annotation (Placement(transformation(
        extent={{-9,-9.5},{9,9.5}},
        rotation=0,
        origin={51,-59.5})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature prescribedTemperature[tube.N_cv] annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=0,
        origin={-10,8})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.ThinWall_L4 thinWall(
    diameter_o=0.55,
    diameter_i=0.5,
    length=tube.length,
    Delta_x=tube.Delta_x,
    N_ax=tube.N_cv,
    T_start=320*ones(tube.N_cv),
    initChoice=ClaRa.Basics.Choices.Init.noInit,
    stateLocation=2) annotation (Placement(transformation(extent={{0,-17},{12,-3}})));
equation
  connect(multiSum.y, massFlowSource.m_flow) annotation (Line(
      points={{69.98,-24},{64,-23},{62,-23}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSink.p, inlet_pressure.y) annotation (Line(
      points={{-56,-35},{-71,-35}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiSum.u[1], mass_flow_1.y) annotation (Line(
      points={{83,-26.1},{83,-24},{88,-24},{88,38},{61,38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(tube.inlet, massFlowSource.steam_a) annotation (Line(
      points={{20,-30},{28,-30},{28,-29},{40,-29}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(massFlowSink.steam_a, tube.outlet) annotation (Line(
      points={{-36,-29},{-10,-29},{-10,-30}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(T_wall.y, realInputMultiplyer.Signal) annotation (Line(
      points={{-69.95,7.5},{-46.42,7.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(mass_flow_2.y, multiSum.u[2]) annotation (Line(
      points={{61,6},{88,6},{88,-24},{83,-24},{83,-21.9}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource.h, inlet_pressure1.y) annotation (Line(
      points={{62,-29},{66,-29},{66,-59.5},{60.9,-59.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realInputMultiplyer.y, prescribedTemperature.T) annotation (Line(
      points={{-33.4,7.5},{-33.4,8},{-17.2,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(prescribedTemperature.port, thinWall.outerPhase) annotation (Line(
      points={{-4,8},{6,8},{6,-3}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(thinWall.innerPhase, tube.heat) annotation (Line(
      points={{6,-17},{6,-25.2},{5,-25.2}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-120},{100,120}}),
                    graphics={Text(
          extent={{-98,112},{100,72}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________
PURPOSE:
test the L4 advanced pipe in a flow reversal scenario to evaluate the numerical robustness and to check for
 physically meaningful behaviour
______________________________________________________________________________________________
",        fontSize=10),Text(
          extent={{-100,120},{100,100}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- 2013-03-05 //JB"),
                      Text(
          extent={{-98,88},{100,48}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________
Scenario:  1) pressure step at inlet (at t=100s 1e5 Pa --> 1.1e5 Pa 
                2) increase of inlet temperature (at t=200s 200e3 J/kg.K --> 220e3 J/kg.K
                3) increase of mass flow (at t = 500s 200kg/s --> 210 kg/s)
                4) increase of outer wall temperature (at t=1000s 293.15K --> 323.15 K)
                5) reversal of mass flow (at t=1500s 210 kg/s --> -100 kg/s)
______________________________________________________________________________________________
",        fontSize=10)}),
    experiment(
      StopTime=2000,
      __Dymola_NumberOfIntervals=1000,
      Tolerance=1e-006,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput(equidistant=false),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=
            true)));
end Test_Pipe_L4_Advanced;
