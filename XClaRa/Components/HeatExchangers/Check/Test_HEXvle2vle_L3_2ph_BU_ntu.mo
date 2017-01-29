within Exergy.XClaRa.Components.HeatExchangers.Check;
model Test_HEXvle2vle_L3_2ph_BU_ntu
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

  HEXvle2vle_L3_2ph_BU_ntu hex(
    redeclare model WallMaterial =
        TILMedia.SolidTypes.TILMedia_Aluminum,
    redeclare model PressureLossTubes =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.VLE_PL.PressureLossCoeffcient_L2
        (                                                                                                    Delta_p_smooth=100, zeta_TOT=5),
    redeclare model PressureLossShell =
        ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3,
    gain_eff=1,
    initTypeTubes=ClaRa.Basics.Choices.Init.noInit,
    m_flow_nom_shell=78,
    p_start_shell=0.023e5,
    CF_geo=1,
    m_flow_nom_tubes=11500,
    p_nom_tubes=1e5,
    h_nom_tubes=60e3,
    h_start_tubes=60e3,
    p_start_tubes=1e5,
    mass_struc=500,
    width_hotwell=2,
    length_hotwell=5,
    level_rel_start=0.2,
    initTypeShell=ClaRa.Basics.Choices.Init.steadyDensity,
    redeclare model HeatTransfer_Shell =
        Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.Constant_L3_ypsDependent
        (                                                                                                    alpha_nom={1000,5000}),
    z_in_tubes=hex.height/2,
    z_out_tubes=hex.height/2,
    z_out_shell=0.05,
    z_in_shell=3.9,
    z_in_aux1=3.9,
    z_in_aux2=3.9,
    redeclare function HeatCapacityAveraging =
        Basics.ControlVolumes.SolidVolumes.Fundamentals.Functions.ArithmeticMean,
    redeclare model HeatTransferTubes =
        Basics.ControlVolumes.Fundamentals.HeatTransport.VLE_HT.NusseltPipe1ph_L2
        (                                                                                                    CF_alpha_tubes=0.5),
    smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments,
    levelOutput=true)
                   annotation (Placement(transformation(extent={{16,-68},{36,-48}})));

  Sensors.Temperature                  Temp_Tubes_out
    annotation (Placement(transformation(extent={{30,-24},{50,-4}})));
  Modelica.Blocks.Sources.Ramp h_steam(
    height=124e3,
    duration=600,
    offset=2212.6e3,
    startTime=10000) annotation (Placement(transformation(extent={{124,32},{104,52}})));
  Modelica.Blocks.Sources.Ramp m_cool(
    duration=100,
    startTime=1000,
    height=0,
    offset=11500) annotation (Placement(transformation(extent={{120,-66},{100,-46}})));
  Modelica.Blocks.Sources.Ramp m_steam(
    startTime=10000,
    offset=76.8,
    duration=600,
    height=-30)   annotation (Placement(transformation(extent={{126,-6},{106,14}})));
  VolumesValvesFittings.Valves.ValveVLE_L1                      valve_shell1(
    checkValve=true,
    redeclare model PressureLoss =
        VolumesValvesFittings.Valves.Fundamentals.QuadraticNominalPoint (                           m_flow_nom=10, Delta_p_nom=250),
    openingInputIsActive=false)
    annotation (Placement(transformation(extent={{-30,-92},{-50,-80}})));
  VolumesValvesFittings.Valves.ValveVLE_L1                      valve_tubes1(
    openingInputIsActive=false,
    checkValve=true,
    redeclare model PressureLoss =
        VolumesValvesFittings.Valves.Fundamentals.LinearNominalPoint (                           Delta_p_nom=1000, m_flow_nom=11500))
    annotation (Placement(transformation(extent={{10,-6},{-10,6}},
        rotation=180,
        origin={54,-32})));
  BoundaryConditions.BoundaryVLE_phxi pressureSink_ph(h_const=300e3, p_const=21e5,
    variable_p=true)                                                               annotation (Placement(transformation(extent={{-74,-96},{-54,-76}})));
  BoundaryConditions.BoundaryVLE_phxi pressureSink_ph1(h_const=2000e3, p_const=250e5,
    variable_p=true)                                                                  annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={84,-32})));
  BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_h(variable_m_flow=true, variable_h=true,
    showData=true)                                                                                 annotation (Placement(transformation(extent={{94,10},{74,-10}})));
  BoundaryConditions.BoundaryVLE_Txim_flow massFlowSource_h1(
    variable_m_flow=true,
    showData=true,
    variable_T=true)
                   annotation (Placement(transformation(extent={{94,-72},{74,-52}})));
  Modelica.Blocks.Sources.Ramp T_cool(
    duration=600,
    offset=13.7 + 273.15,
    startTime=10000,
    height=2)        annotation (Placement(transformation(extent={{120,-96},{100,-76}})));
  inner ClaRa.SimCenter simCenter(
    useHomotopy=true,
    redeclare TILMedia.VLEFluidTypes.TILMedia_SplineWater fluid1,
    showExpertSummary=true)
    annotation (Placement(transformation(extent={{54,30},{74,50}})));
  ClaRa.Visualisation.Hexdisplay_3 hexdisplay_3_1(
    Unit="HEX Temperature in °C",
    y_min=0,
    y_max=100,
    T_o=hex.wall.summary.T_o - fill(273.15, 6),
    T_i=hex.wall.summary.T_i - fill(273.15, 6),
    z_o=hex.wall.summary.eCom.z_o,
    z_i=hex.wall.summary.eCom.z_i)
    annotation (Placement(transformation(extent={{-86,-18},{8,70}})));
  ClaRa.Visualisation.Quadruple quadruple(largeFonts=false)
    annotation (Placement(transformation(extent={{42,-50},{72,-40}})));
  Modelica.Blocks.Sources.Ramp p_steam(
    duration=600,
    startTime=10000,
    height=0.005e5,
    offset=0.023e5) annotation (Placement(transformation(extent={{-100,-90},{-80,-70}})));
  Modelica.Blocks.Sources.Ramp p_cool(
    duration=600,
    startTime=10000,
    height=0,
    offset=1e5) annotation (Placement(transformation(extent={{120,-40},{100,-20}})));
  ClaRa.Visualisation.Quadruple quadruple2(largeFonts=false)
    annotation (Placement(transformation(extent={{38,-75},{68,-65}})));
  ClaRa.Visualisation.Quadruple quadruple1(largeFonts=false)
    annotation (Placement(transformation(extent={{-8,-83},{22,-73}})));
  ClaRa.Visualisation.Quadruple quadruple3(largeFonts=false,
      decimalSpaces(p=3))
    annotation (Placement(transformation(extent={{36,3},{66,13}})));
  ClaRa.Visualisation.DynamicBar level_abs1(
    provideConnector=true,
    u_set=0.8,
    u_high=1,
    u_low=0.6,
    u_max=4,
    u=hex.shell.summary.outline.level_abs)
    annotation (Placement(transformation(extent={{14,-68},{4,-48}})));
  Utilities.Blocks.LimPID PI(
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    Tau_d=60,
    k=0.1,
    u_ref=1,
    y_ref=1,
    y_max=1,
    y_min=0,
    sign=-1,
    y_start=0.5,
    Tau_i=120) annotation (Placement(transformation(extent={{-20,-46},{-30,-36}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=0.8) annotation (Placement(transformation(extent={{2,-46},{-14,-36}})));
  Regression regression(V_liq = hex.shell.summary.outline.volume[1],
    T_shell_out = hex.shell.summary.outlet[1].T,
    p_shell_out = hex.shell.summary.outlet[1].p,
    Q_flow_tot = hex.wall.nTU.summary.Q_flow_tot) annotation (Placement(transformation(extent={{-100,-60},{-80,-40}})));
equation

  connect(valve_tubes1.inlet,Temp_Tubes_out. port) annotation (Line(
      points={{44,-32},{40,-32},{40,-24}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(m_steam.y, massFlowSource_h.m_flow) annotation (Line(
      points={{105,4},{100,4},{100,-6},{96,-6}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource_h.h, h_steam.y) annotation (Line(
      points={{96,0},{98,0},{98,42},{103,42}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(pressureSink_ph.steam_a,valve_shell1. outlet) annotation (Line(
      points={{-54,-86},{-50,-86}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(pressureSink_ph1.steam_a,valve_tubes1. outlet) annotation (Line(
      points={{74,-32},{64,-32}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(hex.Out2, valve_tubes1.inlet) annotation (Line(
      points={{35.8,-52},{40,-52},{40,-32},{44,-32}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(m_cool.y, massFlowSource_h1.m_flow) annotation (Line(
      points={{99,-56},{96,-56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(hex.In1, massFlowSource_h.steam_a) annotation (Line(
      points={{26,-48.2},{26,0},{74,0}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(hex.eye2, quadruple.eye) annotation (Line(points={{37,-50},{42,-50},{42,-45}},            color={190,190,190}));
  connect(p_steam.y, pressureSink_ph.p) annotation (Line(points={{-79,-80},{-74,-80}},           color={0,0,127}));
  connect(pressureSink_ph1.p, p_cool.y) annotation (Line(points={{94,-38},{94,-38},{94,-30},{99,-30}}, color={0,0,127}));
  connect(T_cool.y, massFlowSource_h1.T) annotation (Line(points={{99,-86},{96,-86},{96,-62}}, color={0,0,127}));
  connect(valve_shell1.inlet, hex.Out1) annotation (Line(
      points={{-30,-86},{26,-86},{26,-68}},
      color={0,131,169},
      thickness=0.5));
  connect(massFlowSource_h1.steam_a, hex.In2) annotation (Line(
      points={{74,-62},{35.8,-62}},
      color={0,131,169},
      thickness=0.5));
  connect(quadruple2.eye, massFlowSource_h1.eye) annotation (Line(points={{38,-70},{58,-70},{74,-70}}, color={190,190,190}));
  connect(quadruple1.eye, hex.eye1) annotation (Line(points={{-8,-78},{30,-78},{30,-69}},        color={190,190,190}));
  connect(massFlowSource_h.eye, quadruple3.eye) annotation (Line(points={{74,8},{36,8}}, color={190,190,190}));
  connect(level_abs1.y, PI.u_m) annotation (Line(points={{3,-68},{-2,-68},{-25.05,-68},{-25.05,-47}},
                                                                                          color={0,0,127}));
  connect(PI.y, valve_shell1.opening_in) annotation (Line(points={{-30.5,-41},{-50,-41},{-50,-77},{-40,-77}},  color={0,0,127}));
  connect(realExpression.y, PI.u_s) annotation (Line(points={{-14.8,-41},{-19,-41}}, color={0,0,127}));
  annotation (Diagram(coordinateSystem(extent={{-100,-100},{140,120}},
          preserveAspectRatio=false,
        initialScale=0.1),            graphics={  Text(
          extent={{-100,116},{136,70}},
          lineColor={115,150,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=11,
          textString="______________________________________________________________________________________________
PURPOSE:
>>check HEXvle2vle_L3_2ph_BU_ntu as a condenser in a load change. Test robustness and
prove steady-state initialisation capabilities. Check controlled behaviour.
______________________________________________________________________________________________"),
                       Text(
          extent={{-100,120},{58,102}},
          lineColor={115,150,0},
          fontSize=31,
          textString="TESTED -- 2016-03-02 //TH"),
        Rectangle(
          extent={{-100,120},{140,-100}},
          lineColor={115,150,0},
          lineThickness=0.5)}),                  Icon(coordinateSystem(initialScale=0.1)),
    experiment(
      StopTime=12000,
      __Dymola_NumberOfIntervals=50000,
      Tolerance=1e-005,
      __Dymola_Algorithm="Dassl"),
    __Dymola_experimentSetupOutput(equidistant=false));
end Test_HEXvle2vle_L3_2ph_BU_ntu;
