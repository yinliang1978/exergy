within Exergy.XClaRa.Components.Adapters.Check;
model TestFluidConverters
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb50;
  Fluid2ClaRa fluid2ClaRa
    annotation (Placement(transformation(extent={{-20,34},{0,14}})));
  Modelica.Fluid.Sources.Boundary_ph boundary(
    use_p_in=false,
    use_h_in=false,
    nPorts=1,
    redeclare package Medium = Modelica.Media.Water.StandardWater,
    p=1000000,
    h=100e3)
    annotation (Placement(transformation(extent={{-144,14},{-124,34}})));
  BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_T(variable_m_flow=true, variable_h=true,
    showData=true)                                                                                 annotation (Placement(transformation(extent={{60,14},{40,34}})));
  Modelica.Blocks.Sources.Ramp ramp(
    height=-2,
    duration=0.1,
    offset=1,
    startTime=0.75) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={90,30})));
  Modelica.Blocks.Sources.Ramp ramp1(
    duration=0.1,
    offset=1e5,
    startTime=0.25,
    height=3e6)     annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={90,0})));
  Modelica.Fluid.Sensors.MassFlowRate massFlowRate(redeclare package Medium =
        Modelica.Media.Water.StandardWater)
    annotation (Placement(transformation(extent={{-98,14},{-118,34}})));
  Modelica.Fluid.Sensors.Temperature temperature(redeclare package Medium =
        Modelica.Media.Water.StandardWater)
    annotation (Placement(transformation(extent={{-30,24},{-50,44}})));
  inner ClaRa.SimCenter simCenter annotation (Placement(
        transformation(extent={{70,-214},{90,-194}})));
  BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_T1(variable_m_flow=true, variable_h=true,
    showData=true)                                                                                  annotation (Placement(transformation(extent={{62,-68},{42,-48}})));
  ClaRa.Visualisation.Quadruple quadruple1 annotation (Placement(
        transformation(extent={{30,-90},{-8,-80}})));
  ThermoPower2ClaRa thermoPower2ClaRa
    annotation (Placement(transformation(extent={{-32,-48},{-12,-68}})));
  ClaRa.Visualisation.DynDisplay dynDisplay2(
    varname="MSL.T",
    unit="°C",
    x1=sensT.T - 273.15)
    annotation (Placement(transformation(extent={{-60,8},{-24,20}})));
  ClaRa.Visualisation.DynDisplay dynDisplay3(
    unit="kg/s",
    x1=sensT.outlet.w,
    varname="ThermoPower.m_flow") annotation (Placement(
        transformation(extent={{-80,-72},{-44,-60}})));
  Fundamentals.SourceP sourceP(p0=1000000) annotation (Placement(transformation(extent={{-144,-68},{-124,-48}})));
  Fundamentals.SensT sensT annotation (Placement(transformation(extent={{-110,-64},{-90,-44}})));
  BoundaryConditions.BoundaryVLE_hxim_flow massFlowSource_T2(variable_m_flow=true, variable_h=true,
    showData=true)                                                                                  annotation (Placement(transformation(extent={{-136,-162},{-116,-142}})));
  ClaRa.Visualisation.Quadruple quadruple2 annotation (Placement(
        transformation(extent={{-110,-165},{-72,-155}})));
  ClaRa2ThermoPower thermoPower2ClaRa1
    annotation (Placement(transformation(extent={{-30,-162},{-10,-142}})));
  Fundamentals.SinkP sourceP1(p0=1000000) annotation (Placement(transformation(extent={{70,-162},{90,-142}})));
  Fundamentals.SensT sensT1 annotation (Placement(transformation(extent={{-8,-158},{12,-138}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue(significantDigits=3)
                                                       annotation (Placement(transformation(extent={{-80,24},{-100,44}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue1(significantDigits=3)
                                                        annotation (Placement(transformation(extent={{-114,34},{-134,54}})));
  Sensors.vleMassflowSensor Flow annotation (Placement(transformation(extent={{34,24},{14,44}})));
  Sensors.Temperature Temperature annotation (Placement(transformation(extent={{24,24},{4,4}})));
  Modelica.Blocks.Math.Feedback feedback annotation (Placement(transformation(extent={{-54,44},{-74,24}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=273.15) annotation (Placement(transformation(extent={{-100,40},{-80,60}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue2(significantDigits=4)
                                                       annotation (Placement(transformation(extent={{-60,-58},{-40,-38}})));
  Modelica.Blocks.Math.Feedback feedback1
                                         annotation (Placement(transformation(extent={{-86,-38},{-66,-58}})));
  Modelica.Blocks.Sources.RealExpression realExpression1(
                                                        y=273.15) annotation (Placement(transformation(extent={{-116,-40},{-96,-20}})));
  Sensors.Temperature Temperature1(unitOption=2)
                                  annotation (Placement(transformation(extent={{10,-58},{-10,-38}})));
  Sensors.vleMassflowSensor Flow1
                                 annotation (Placement(transformation(extent={{34,-58},{14,-38}})));
  Modelica.Blocks.Interaction.Show.RealValue realValue3(significantDigits=4)
                                                       annotation (Placement(transformation(extent={{46,-152},{66,-132}})));
  Modelica.Blocks.Math.Feedback feedback2
                                         annotation (Placement(transformation(extent={{22,-132},{42,-152}})));
  Modelica.Blocks.Sources.RealExpression realExpression2(
                                                        y=273.15) annotation (Placement(transformation(extent={{66,-136},{46,-116}})));
  Sensors.vleMassflowSensor Flow2
                                 annotation (Placement(transformation(extent={{-106,-152},{-86,-132}})));
  Sensors.Temperature Temperature2(unitOption=2)
                                  annotation (Placement(transformation(extent={{-78,-152},{-58,-132}})));
  inner Modelica.Fluid.System system annotation (Placement(transformation(extent={{-150,-214},{-130,-194}})));
equation
  connect(ramp.y, massFlowSource_T.m_flow) annotation (Line(
      points={{79,30},{62,30}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp1.y, massFlowSource_T.h) annotation (Line(
      points={{79,1.33227e-015},{76,1.33227e-015},{76,24},{62,24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowRate.port_a, fluid2ClaRa.port_a) annotation (Line(
      points={{-98,24},{-98,24},{-19.8,24}},
      color={0,127,255},
      smooth=Smooth.Bezier));
  connect(boundary.ports[1], massFlowRate.port_b) annotation (Line(
      points={{-124,24},{-118,24}},
      color={0,127,255},
      smooth=Smooth.Bezier));
  connect(temperature.port, massFlowRate.port_a) annotation (Line(
      points={{-40,24},{-98,24}},
      color={0,127,255},
      smooth=Smooth.Bezier));
  connect(ramp.y, massFlowSource_T1.m_flow) annotation (Line(
      points={{79,30},{74,30},{74,-52},{64,-52}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp1.y, massFlowSource_T1.h) annotation (Line(
      points={{79,1.33227e-015},{76,1.33227e-015},{76,-58},{64,-58}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(massFlowSource_T1.eye, quadruple1.eye) annotation (Line(
      points={{42,-66},{36,-66},{36,-85},{30,-85}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(sensT.outlet, thermoPower2ClaRa.flangeA) annotation (Line(
      points={{-94,-58},{-31.8,-58}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(sourceP.flange, sensT.inlet) annotation (Line(
      points={{-124,-58},{-106,-58}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(massFlowSource_T2.eye,quadruple2. eye) annotation (Line(
      points={{-116,-160},{-110,-160}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(thermoPower2ClaRa1.flangeB, sensT1.inlet) annotation (Line(
      points={{-10,-152},{-4,-152}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(ramp1.y, massFlowSource_T2.h) annotation (Line(
      points={{79,1.33227e-015},{76,1.33227e-015},{76,-108},{-146,-108},{-146,-152},{-138,-152}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(ramp.y, massFlowSource_T2.m_flow) annotation (Line(
      points={{79,30},{74,30},{74,-108},{-146,-108},{-146,-146},{-138,-146}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(realValue1.numberPort, massFlowRate.m_flow) annotation (Line(points={{-112.5,44},{-108,44},{-108,35}},    color={0,0,127}));
  connect(Flow.inlet, massFlowSource_T.steam_a) annotation (Line(
      points={{34,24},{40,24}},
      color={0,131,169},
      thickness=0.5));
  connect(Flow.outlet, fluid2ClaRa.steam_a) annotation (Line(
      points={{14,24},{-0.05,24},{-0.05,24.0025}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(fluid2ClaRa.steam_a, Temperature.port) annotation (Line(
      points={{-0.05,24.0025},{14,24.0025},{14,24}},
      color={0,131,169},
      thickness=0.5));
  connect(temperature.T, feedback.u1) annotation (Line(points={{-47,34},{-47,34},{-56,34}}, color={0,0,127}));
  connect(feedback.y, realValue.numberPort) annotation (Line(points={{-73,34},{-76,34},{-78.5,34}}, color={0,0,127}));
  connect(realExpression.y, feedback.u2) annotation (Line(points={{-79,50},{-64,50},{-64,42}}, color={0,0,127}));
  connect(realExpression1.y, feedback1.u2) annotation (Line(points={{-95,-30},{-76,-30},{-76,-40}}, color={0,0,127}));
  connect(Temperature1.port, Flow1.outlet) annotation (Line(
      points={{0,-58},{14,-58}},
      color={0,131,169},
      thickness=0.5));
  connect(Flow1.outlet, thermoPower2ClaRa.outlet) annotation (Line(
      points={{14,-58},{-12,-58}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(Flow1.inlet, massFlowSource_T1.steam_a) annotation (Line(
      points={{34,-58},{42,-58}},
      color={0,131,169},
      thickness=0.5));
  connect(sensT.T, feedback1.u1) annotation (Line(points={{-92,-48},{-88,-48},{-84,-48}}, color={0,0,127}));
  connect(feedback1.y, realValue2.numberPort) annotation (Line(points={{-67,-48},{-61.5,-48},{-61.5,-48}}, color={0,0,127}));
  connect(feedback2.u1, sensT1.T) annotation (Line(points={{24,-142},{10,-142}}, color={0,0,127}));
  connect(feedback2.y, realValue3.numberPort) annotation (Line(points={{41,-142},{44.5,-142}}, color={0,0,127}));
  connect(realExpression2.y, feedback2.u2) annotation (Line(points={{45,-126},{32,-126},{32,-134}}, color={0,0,127}));
  connect(massFlowSource_T2.steam_a, Flow2.inlet) annotation (Line(
      points={{-116,-152},{-106,-152}},
      color={0,131,169},
      thickness=0.5));
  connect(Flow2.outlet, Temperature2.port) annotation (Line(
      points={{-86,-152},{-86,-152},{-68,-152}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5));
  connect(Temperature2.port, thermoPower2ClaRa1.inlet) annotation (Line(
      points={{-68,-152},{-38,-152},{-30,-152}},
      color={0,131,169},
      thickness=0.5));
  connect(sensT1.outlet, sourceP1.flange) annotation (Line(points={{8,-152},{70,-152},{70,-152}}, color={0,0,255}));
  annotation (Diagram(coordinateSystem(extent={{-150,-200},{100,100}}, initialScale=0.1), graphics={
        Polygon(
          points={{-150,80},{36,80},{-28,0},{-150,0},{-150,80}},
          lineColor={215,215,215},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{36,80},{100,80},{100,0},{-28,0},{36,80}},
          lineColor={215,215,215},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-32,78},{28,70}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="Modelica Standard Library"),
        Text(
          extent={{20,78},{80,70}},
          lineColor={0,131,169},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="ClaRa Library",
          fontName="Miso"),
        Polygon(
          points={{-150,-16},{20,-16},{-60,-96},{-150,-96},{-150,-16}},
          lineColor={215,215,215},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{20,-16},{100,-16},{100,-96},{-60,-96},{20,-16}},
          lineColor={215,215,215},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{4,-18},{64,-26}},
          lineColor={0,131,169},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="ClaRa Library",
          fontName="Miso"),
        Text(
          extent={{-40,-18},{12,-26}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="ThermoPower Library"),
        Polygon(
          points={{-52,-112},{100,-112},{100,-192},{-60,-192},{-52,-112}},
          lineColor={215,215,215},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-152,-112},{14,-112},{-54,-192},{-152,-192},{-152,-112}},
          lineColor={215,215,215},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-46,-182},{6,-190}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="ThermoPower Library"),
        Text(
          extent={{-102,-180},{-42,-188}},
          lineColor={0,131,169},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="ClaRa Library",
          fontName="Miso")}), Icon(coordinateSystem(                                initialScale=0.1)));
end TestFluidConverters;
