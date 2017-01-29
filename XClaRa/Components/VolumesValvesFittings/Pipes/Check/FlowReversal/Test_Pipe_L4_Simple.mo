within Exergy.XClaRa.Components.VolumesValvesFittings.Pipes.Check.FlowReversal;
model Test_Pipe_L4_Simple
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
        origin={77,-38})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource(
    m_flow_const=0.1,
    variable_m_flow=true,
    h_const=200e3,
    m_flow_nom=0,
    variable_h=true,
    p_nom=1000) annotation (Placement(transformation(extent={{58,
            -53},{38,-33}})));
  inner ClaRa.SimCenter simCenter(redeclare replaceable
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1,
      useHomotopy=true) annotation (Placement(transformation(
          extent={{-80,-114},{-60,-94}})));
  PipeFlowVLE_L4_Simple tube(
    z_in=0,
    z_out=0,
    showExpertSummary=true,
    showData=true,
    m_flow_nom=100,
    diameter_i=0.5,
    h_start=ones(tube.N_cv)*200e3,
    initType=ClaRa.Basics.Choices.Init.steadyState,
    length=50,
    redeclare model PressureLoss =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L4,
    p_start=linspace(
        3.2e5,
        3.0e5,
        tube.N_cv),
    N_cv=50,
    Delta_x=ones(tube.N_cv)*tube.length/tube.N_cv,
    frictionAtInlet=true,
    frictionAtOutlet=true) annotation (Placement(transformation(extent={{20,-49},{-16,-36}})));

  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_phxi massFlowSink(
    variable_p=true,
    m_flow_nom=100,
    p_const=1000000,
    h_const=200e3,
    Delta_p=100000) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-48,-43})));
  inner Modelica.Fluid.System system
    annotation (Placement(transformation(extent={{60,-114},{80,-94}})));
  Modelica.Blocks.Sources.Step inlet_pressure(
    offset=1e5,
    startTime=100,
    height=1e4)
    annotation (Placement(transformation(extent={{-92,-59},{-72,-39}})));
  Modelica.Blocks.Sources.Ramp mass_flow_1(
    duration=1,
    height=10,
    offset=100,
    startTime=500) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={48,22})));

  Modelica.Blocks.Sources.Ramp T_wall(
    offset=293.15,
    startTime=1000,
    duration=100,
    height=30) annotation (Placement(transformation(
        extent={{-10.5,-10.5},{10.5,10.5}},
        rotation=0,
        origin={-83.5,-0.5})));
  Utilities.Blocks.RealInputMultiplyer realInputMultiplyer(N=tube.N_cv) annotation (Placement(transformation(extent={{-60,-10},{-46,9}})));

  Modelica.Blocks.Sources.Ramp mass_flow_2(
    offset=100,
    startTime=1500,
    height=-300,
    duration=200) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={48,-10})));

  Modelica.Blocks.Sources.Step inlet_pressure1(
    height=20e3,
    offset=200e3,
    startTime=200) annotation (Placement(transformation(
        extent={{-9,-9.5},{9,9.5}},
        rotation=0,
        origin={49,-71.5})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature prescribedTemperature[tube.N_cv] annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=0,
        origin={-22,0})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.ThinWall_L4 thinWall(
    diameter_o=0.55,
    diameter_i=0.5,
    length=tube.length,
    Delta_x=tube.Delta_x,
    initChoice=ClaRa.Basics.Choices.Init.noInit,
    stateLocation=2,
    N_ax=tube.N_cv,
    T_start=320*ones(tube.N_cv)) annotation (Placement(transformation(extent={{-4,-23},{8,-9}})));
equation
  connect(multiSum.y, massFlowSource.m_flow) annotation (Line(
      points={{69.98,-38},{60,-38},{60,-37}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSink.p, inlet_pressure.y) annotation (Line(
      points={{-58,-49},{-71,-49}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(tube.inlet, massFlowSource.steam_a) annotation (Line(
      points={{20,-42.5},{26,-42.5},{26,-43},{38,-43}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(massFlowSink.steam_a, tube.outlet) annotation (Line(
      points={{-38,-43},{-24,-43},{-24,-42.5},{-16,-42.5}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(T_wall.y, realInputMultiplyer.Signal) annotation (Line(
      points={{-71.95,-0.5},{-60.42,-0.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource.h, inlet_pressure1.y) annotation (Line(
      points={{60,-43},{66,-43},{66,-71.5},{58.9,-71.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realInputMultiplyer.y, prescribedTemperature.T) annotation (Line(
      points={{-47.4,-0.5},{-47.4,0},{-29.2,0}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(prescribedTemperature.port, thinWall.outerPhase) annotation (Line(
      points={{-16,0},{2,0},{2,-9}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(thinWall.innerPhase, tube.heat) annotation (Line(
      points={{2,-23},{2,-37.3}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(mass_flow_1.y, multiSum.u[1]) annotation (Line(points={{59,22},{94,22},{94,-40.1},{83,-40.1}}, color={0,0,127}));
  connect(mass_flow_2.y, multiSum.u[2]) annotation (Line(points={{59,-10},{88,-10},{88,-35.9},{83,-35.9}}, color={0,0,127}));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-120},{100,120}}),
                    graphics={Text(
          extent={{-98,112},{100,72}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________
PURPOSE:
test the L4 simple pipe in a flow reversal scenario to evaluate the numerical robustness and to check for
 physically meaningful behaviour
______________________________________________________________________________________________
",        fontSize=10),Text(
          extent={{-100,120},{100,100}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- 2013-03-05 //JB"),Text(
          extent={{-98,66},{80,24}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________________________
Remarks: - tube is initialised steady state in pressure and enthalpy
______________________________________________________________________________________________________________
",        fontSize=8),Text(
          extent={{-98,86},{100,46}},
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
      StopTime=3000,
      __Dymola_NumberOfIntervals=1000,
      Tolerance=1e-006,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput,
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=
            true)));
end Test_Pipe_L4_Simple;
