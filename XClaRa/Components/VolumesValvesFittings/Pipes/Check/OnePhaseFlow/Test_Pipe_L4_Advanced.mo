within Exergy.XClaRa.Components.VolumesValvesFittings.Pipes.Check.OnePhaseFlow;
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
        origin={73,-54})));
  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_hxim_flow
    massFlowSource(
    m_flow_const=0.1,
    variable_m_flow=true,
    h_const=200e3,
    m_flow_nom=0,
    variable_h=true,
    p_nom=1000) annotation (Placement(transformation(extent={{60,
            -69},{40,-49}})));
  inner ClaRa.SimCenter simCenter(redeclare replaceable
      TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater fluid1,
      useHomotopy=true) annotation (Placement(transformation(
          extent={{-100,-140},{-80,-120}})));
  PipeFlowVLE_L4_Advanced tube(
    z_in=0,
    z_out=0,
    showExpertSummary=true,
    showData=true,
    m_flow_nom=100,
    diameter_i=0.5,
    Delta_x=ones(tube.N_cv)*tube.length/tube.N_cv,
    h_start=ones(tube.N_cv)*200e3,
    length=50,
    N_cv=50,
    redeclare model PressureLoss =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearPressureLoss_L4,
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L4,
    p_start=ones(tube.N_cv)*1e5,
    frictionAtInlet=true,
    frictionAtOutlet=true,
    initType=ClaRa.Basics.Choices.Init.steadyState,
    suppressHighFrequencyOscillations=true) annotation (Placement(transformation(extent={{14,-63},{-8,-55}})));

  Exergy.XClaRa.Components.BoundaryConditions.BoundaryVLE_phxi massFlowSink(
    variable_p=true,
    m_flow_nom=100,
    p_const=1000000,
    Delta_p=100000,
    h_const=200e3) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=180,
        origin={-46,-59})));
  inner Modelica.Fluid.System system
    annotation (Placement(transformation(extent={{80,-140},{100,-120}})));
  Modelica.Blocks.Sources.Step inlet_pressure(
    offset=1e5,
    startTime=100,
    height=1e4)
    annotation (Placement(transformation(extent={{-100,-75},{-80,-55}})));
  Modelica.Blocks.Sources.Ramp mass_flow_1(
    duration=1,
    height=10,
    offset=100,
    startTime=500) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={50,-10})));

  Modelica.Blocks.Sources.Ramp T_wall(
    duration=4,
    startTime=1000,
    offset=293.15,
    height=200) annotation (Placement(transformation(
        extent={{-10.5,-10.5},{10.5,10.5}},
        rotation=0,
        origin={-89.5,-22.5})));
  Utilities.Blocks.RealInputMultiplyer realInputMultiplyer(N=tube.N_cv) annotation (Placement(transformation(extent={{-64,-32},{-50,-13}})));

  Modelica.Blocks.Sources.Ramp mass_flow_2(
    duration=1,
    offset=100,
    height=-50,
    startTime=1500) annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=0,
        origin={90,-10})));

  Modelica.Blocks.Sources.Step inlet_pressure1(
    height=20e3,
    offset=200e3,
    startTime=2000) annotation (Placement(transformation(
        extent={{9,-9.5},{-9,9.5}},
        rotation=0,
        origin={91,-89.5})));

  Modelica.Thermal.HeatTransfer.Sources.PrescribedTemperature prescribedTemperature[tube.N_cv] annotation (Placement(transformation(
        extent={{-6,-6},{6,6}},
        rotation=0,
        origin={-34,-22})));
  ClaRa.Basics.ControlVolumes.SolidVolumes.ThinWall_L4 thinWall1(
    diameter_o=0.55,
    diameter_i=0.5,
    length=tube.length,
    Delta_x=tube.Delta_x,
    N_ax=tube.N_cv,
    T_start=320*ones(tube.N_cv),
    initChoice=ClaRa.Basics.Choices.Init.noInit,
    stateLocation=2) annotation (Placement(transformation(extent=
            {{-10,-42},{16,-32}})));
equation
  connect(multiSum.y, massFlowSource.m_flow) annotation (Line(
      points={{65.98,-54},{64,-54},{62,-53}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSink.p, inlet_pressure.y) annotation (Line(
      points={{-56,-65},{-79,-65}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(multiSum.u[1], mass_flow_1.y) annotation (Line(
      points={{79,-56.1},{79,-32},{69,-32},{69,-10},{61,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(tube.inlet, massFlowSource.steam_a) annotation (Line(
      points={{14,-59},{40,-59}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(massFlowSink.steam_a, tube.outlet) annotation (Line(
      points={{-36,-59},{-8,-59}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(T_wall.y, realInputMultiplyer.Signal) annotation (Line(
      points={{-77.95,-22.5},{-64.42,-22.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(mass_flow_2.y, multiSum.u[2]) annotation (Line(
      points={{79,-10},{79,-51.9}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource.h, inlet_pressure1.y) annotation (Line(
      points={{62,-59},{66,-59},{66,-89.5},{81.1,-89.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realInputMultiplyer.y, prescribedTemperature.T) annotation (Line(
      points={{-51.4,-22.5},{-51.4,-22},{-41.2,-22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(thinWall1.innerPhase, tube.heat) annotation (Line(
      points={{3,-42},{3,-55.8}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(thinWall1.outerPhase, prescribedTemperature.port) annotation (Line(
      points={{3,-32},{3,-22},{-28,-22}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-140},{100,
            120}}), graphics={Text(
          extent={{-100,112},{98,72}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________
PURPOSE:
test the advanced L4 pipe in a various number of steps concerning mass flow rate heat flow rate 
and inlet enthalpy to evaluate the numerical robustness and to check for physically meaningful behaviour
______________________________________________________________________________________________
",        fontSize=10),Text(
          extent={{-100,120},{100,100}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- 2013-04-18 //JB"),Text(
          extent={{-98,58},{80,16}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=7,
          textString="______________________________________________________________________________________________________________
Remarks: - try the effect of the supressHighFrequencyOscillations at |Physical Effects|Momentum Balance
                                          --> set this to false in order to get agreement to Modelica Fluid pipe tester
        
        
______________________________________________________________________________________________________________"),
          Text(
          extent={{-100,88},{98,48}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________
Scenario:  1) pressure step at inlet (at t=100s 1e5 Pa --> 1.1e5 Pa 
                2) increase of mass flow (at t = 500s 200kg/s --> 210 kg/s)
                3) increase of outer wall temperature (at t=1000s 293.15K --> 493.15 K)
                4) decrease of mass flow (at t=1500s 210 kg/s --> 160 kg/s)
                5) increase of inlet temperature (at t=2000s 200e3 J/kg.K --> 220e3 J/kg.K
______________________________________________________________________________________________
",        fontSize=10)}),
    experiment(
      StopTime=3000,
      __Dymola_NumberOfIntervals=10000,
      Tolerance=1e-006,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput(equdistant=false, events=false),
    Icon(coordinateSystem(extent={{-100,-100},{100,100}}, preserveAspectRatio=
            true)));
end Test_Pipe_L4_Advanced;
