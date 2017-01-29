within Exergy.XClaRa.Components.HeatExchangers.Check;
model Test_HEXvle2vle_L3_2ph_CH_simple
 extends ClaRa.Basics.Icons.PackageIcons.ExecutableRegressiong100;
model Regression
  extends ClaRa.Basics.Icons.RegressionSummary;
  Modelica.Blocks.Interfaces.RealInput V_liq "Liquid shell volume";
  Modelica.Blocks.Interfaces.RealInput T_shell_out "Shell outlet temperature";
  Modelica.Blocks.Interfaces.RealInput p_shell_out "IP turbine outlet enthalpy";
  Modelica.Blocks.Interfaces.RealInput Q_flow_tot "Total heat flow";

  Real y_Q_flow_tot_int = integrator1.y;
  Real y_Q_flow_tot = Q_flow_tot;

  Real y_V_liq_int = integrator2.y;
  Real y_V_liq = V_liq;

  Real y_T_shell_out_int = integrator3.y;
  Real y_T_shell_out = T_shell_out;
  Real y_p_shell_out_int = integrator4.y;
  Real y_p_shell_out = p_shell_out;

  protected
  Utilities.Blocks.Integrator integrator1(u=Q_flow_tot, startTime=
          1000);
  Utilities.Blocks.Integrator integrator2(u=V_liq, startTime=1000);
  Utilities.Blocks.Integrator integrator3(u=T_shell_out, startTime=
          1000);
  Utilities.Blocks.Integrator integrator4(u=p_shell_out, startTime=
          1000);
end Regression;

  HEXvle2vle_L3_2ph_CH_simple hex(
    mass_struc=1,
    redeclare model WallMaterial =
        TILMedia.SolidTypes.TILMedia_Aluminum,
    redeclare model PressureLossTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.QuadraticNominalPoint_L2,
    redeclare model PressureLossShell =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.QuadraticParallelZones_L3,
    z_in_shell=10,
    redeclare model HeatTransferTubes =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L2
        (                                                                                                    alpha_nom=5000),
    T_w_start=linspace(
        200 + 273.15,
        450 + 273.15,
        3),
    p_start_tubes=250e5,
    initTypeWall=ClaRa.Basics.Choices.Init.steadyState,
    m_flow_nom_shell=42,
    p_nom_shell=53e5,
    h_nom_shell=3000e3,
    p_start_shell=53e5,
    N_passes=3,
    z_in_tubes=0.1,
    z_out_tubes=10,
    p_nom_tubes=250e5,
    h_nom_tubes=1000e3,
    h_start_tubes=1000e3,
    m_flow_nom_tubes=416,
    z_out_shell=0.1,
    initTypeTubes=ClaRa.Basics.Choices.Init.noInit,
    redeclare model HeatTransfer_Shell =
        Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.Constant_L3_ypsDependent
        (                                                                                                    alpha_nom={1000,3000}),
    level_rel_start=0.2,
    N_tubes=500)                                                                                                     annotation (Placement(transformation(extent={{-6,-72},{14,-52}})));

  Sensors.Temperature                  Temp_Shell_in
    annotation (Placement(transformation(extent={{14,-22},{34,-42}})));
  Sensors.Temperature                  Temp_Tubes_in
    annotation (Placement(transformation(extent={{38,-68},{58,-88}})));
  Modelica.Blocks.Sources.Ramp h_hot(
    offset=2942e3,
    duration=600,
    startTime=10000,
    height=80e3) annotation (Placement(transformation(extent={{100,-32},{80,-12}})));
  Modelica.Blocks.Sources.Ramp m_cold(
    duration=600,
    offset=416,
    startTime=10000,
    height=-166) annotation (Placement(transformation(extent={{100,-64},{80,-44}})));
  Modelica.Blocks.Sources.Ramp m_hot(
    startTime=10000,
    duration=600,
    offset=42.7,
    height=-20) annotation (Placement(transformation(extent={{100,0},{80,20}})));
  VolumesValvesFittings.Valves.ValveVLE_L1                      valve_shell1(
    checkValve=true,
    redeclare model PressureLoss =
        VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (                           m_flow_nom=50, Delta_p_nom=2e5),
    openingInputIsActive=false)
    annotation (Placement(transformation(extent={{-48,-94},{-68,-82}})));
  BoundaryConditions.BoundaryVLE_phxi pressureSink_ph(h_const=300e3, p_const=2100000,
    variable_p=true)                                                                  annotation (Placement(transformation(extent={{-92,-98},{-72,-78}})));
  BoundaryConditions.BoundaryVLE_phxi pressureSink_ph1(h_const=2000e3, p_const=25000000,
    variable_p=true)                                                                     annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-82,-60})));
  BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_h(variable_m_flow=true, variable_h=true) annotation (Placement(transformation(extent={{56,-32},{36,-12}})));
  BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_h1(variable_m_flow=true, variable_h=true) annotation (Placement(transformation(extent={{68,-70},{48,-50}})));
  Modelica.Blocks.Sources.Ramp h_cold(
    duration=600,
    offset=961e3,
    startTime=10000,
    height=-106e3) annotation (Placement(transformation(extent={{100,-94},{80,-74}})));
  inner ClaRa.SimCenter simCenter(
    useHomotopy=true,
    redeclare TILMedia.VLEFluidTypes.TILMedia_InterpolatedWater
      fluid1,
    showExpertSummary=true)
    annotation (Placement(transformation(extent={{40,40},{80,60}})));
  ClaRa.Visualisation.Hexdisplay_3 hexdisplay_3_1(
    T_o={hex.shell.summary.inlet[1].T,hex.shell.summary.outlet[1].T,
        hex.shell.summary.outlet[1].T,hex.shell.summary.outlet[1].T,
        hex.shell.summary.outlet[1].T,hex.shell.summary.outlet[1].T},
    T_i={hex.tubes.summary.inlet.T,hex.tubes.summary.outlet.T,hex.tubes.summary.outlet.T,
        hex.tubes.summary.outlet.T,hex.tubes.summary.outlet.T,hex.tubes.summary.outlet.T},
    Unit="HEX Temperature in K",
    y_min=500,
    y_max=800,
    z_i={0,1,1,1,1,1},
    z_o={0,1,1,1,1,1}) annotation (Placement(transformation(extent={{
            -120,-32},{-26,56}})));

  ClaRa.Visualisation.Quadruple quadruple(largeFonts=false)
    annotation (Placement(transformation(extent={{-42,-50},{-10,-34}})));
  ClaRa.Visualisation.Quadruple quadruple1(largeFonts=false)
    annotation (Placement(transformation(extent={{12,-100},{44,-84}})));
  Modelica.Blocks.Sources.Ramp p_cold(
    duration=600,
    startTime=10000,
    offset=250e5,
    height=-95e5) annotation (Placement(transformation(extent={{-120,-64},{-100,-44}})));
  Modelica.Blocks.Sources.Ramp p_hot(
    duration=600,
    startTime=10000,
    offset=53e5,
    height=-20e5) annotation (Placement(transformation(extent={{-120,-92},{-100,-72}})));
  ClaRa.Visualisation.DynamicBar level_abs1(
    provideConnector=true,
    u=hex.shell.summary.outline.level_abs,
    u_set=2,
    u_high=3,
    u_low=1,
    u_max=10) annotation (Placement(transformation(extent={{-8,-72},{
            -18,-52}})));
  Utilities.Blocks.LimPID PI(
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    Tau_d=60,
    k=0.1,
    u_ref=1,
    y_ref=1,
    y_max=1,
    y_min=0,
    y_start=0.5,
    Tau_i=120,
    sign=1) annotation (Placement(transformation(extent={{-28,-77},{-38,-67}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=2) annotation (Placement(transformation(extent={{-8,-86},{-24,-76}})));

     Regression regression(V_liq = hex.shell.summary.outline.volume[1],
     T_shell_out = hex.shell.summary.outlet[1].T,
     p_shell_out = hex.shell.summary.outlet[1].p,
     Q_flow_tot = -hex.summary.outline.Q_flow) annotation (Placement(transformation(extent={{-120,80},{-100,100}})));
equation

  connect(m_hot.y, massFlowSource_h.m_flow) annotation (Line(
      points={{79,10},{64,10},{64,-16},{58,-16}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource_h.h, h_hot.y) annotation (Line(
      points={{58,-22},{79,-22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pressureSink_ph.steam_a,valve_shell1. outlet) annotation (Line(
      points={{-72,-88},{-68,-88}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(massFlowSource_h1.h, h_cold.y) annotation (Line(
      points={{70,-60},{70,-84},{79,-84}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(m_cold.y, massFlowSource_h1.m_flow) annotation (Line(
      points={{79,-54},{70,-54}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(hex.In1, massFlowSource_h.steam_a) annotation (Line(
      points={{4,-52.2},{4,-22},{36,-22}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(hex.In2, massFlowSource_h1.steam_a) annotation (Line(
      points={{14,-60},{48,-60}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(hex.In1, Temp_Shell_in.port) annotation (Line(
      points={{4,-52.2},{4,-22},{24,-22}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(hex.In2, Temp_Tubes_in.port) annotation (Line(
      points={{14,-60},{49,-60},{49,-68},{48,-68}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(hex.eye1, quadruple1.eye) annotation (Line(points={{8,-73},{8,-92},{12,-92}},       color={190,190,190}));
  connect(hex.eye2, quadruple.eye) annotation (Line(points={{-7,-62},{-8,-62},{-8,-42},{-42,-42}},   color={190,190,190}));
  connect(p_cold.y, pressureSink_ph1.p) annotation (Line(points={{-99,-54},{-96,-54},{-92,-54}}, color={0,0,127}));
  connect(p_hot.y, pressureSink_ph.p) annotation (Line(points={{-99,-82},{-92,-82}}, color={0,0,127}));
  connect(level_abs1.y, PI.u_s) annotation (Line(points={{-19,-72},{-22,-72},{-27,-72}}, color={0,0,127}));
  connect(realExpression.y, PI.u_m) annotation (Line(points={{-24.8,-81},{-33.05,-81},{-33.05,-78}},
                                                                                          color={0,0,127}));
  connect(hex.Out1, valve_shell1.inlet) annotation (Line(
      points={{4,-72},{4,-88},{-48,-88}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(PI.y, valve_shell1.opening_in) annotation (Line(points={{-38.5,-72},{-58,-72},{-58,-79}},  color={0,0,127}));
  connect(hex.Out2, pressureSink_ph1.steam_a) annotation (Line(
      points={{-6,-60},{-72,-60}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  annotation (Diagram(coordinateSystem(extent={{-120,-100},{100,100}},
          preserveAspectRatio=false,
        initialScale=0.1),            graphics={  Text(
          extent={{-96,102},{142,56}},
          lineColor={115,150,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=11,
          textString="______________________________________________________________________________________________
PURPOSE:
>>check HEXvle2vle_L3_2ph_CH_simple as a high pressure preheater in a load change. 
Test robustness and prove steady-state initialisation capabilities. Check controlled and uncontrolled behaviour.
______________________________________________________________________________________________"),
                       Text(
          extent={{-114,102},{44,84}},
          lineColor={115,150,0},
          fontSize=30,
          textString="TESTED -- 2014-10-16 //TH"),
        Rectangle(
          extent={{-120,100},{100,-100}},
          lineColor={115,150,0},
          lineThickness=0.5)}),                  Icon(coordinateSystem(extent={{-100,-100},{100,100}})),
    experiment(StopTime=12000, Tolerance=1e-005),
    __Dymola_experimentSetupOutput);
end Test_HEXvle2vle_L3_2ph_CH_simple;
