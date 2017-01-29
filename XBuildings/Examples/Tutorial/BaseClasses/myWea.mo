within Exergy.XBuildings.Examples.Tutorial.BaseClasses;
model myWea

  Buildings.BoundaryConditions.WeatherData.Bus
      weaBus "Weather data bus" annotation (Placement(transformation(extent={{58,-12},
            {78,8}}),            iconTransformation(extent={{58,-12},
            {78,8}})));
  Modelica.Blocks.Sources.RealExpression tem(y=273.15 + 33 + 8*sin(2*
        3.14*time/(24*3600)))
    annotation (Placement(transformation(extent={{-50,50},{-30,70}})));
  Modelica.Blocks.Sources.RealExpression pressure(y=101325)
    annotation (Placement(transformation(extent={{-48,-38},{-28,-18}})));
  Modelica.Blocks.Sources.RealExpression relativeHum(y=0.9)
    annotation (Placement(transformation(extent={{-44,-66},{-24,-46}})));
equation
  connect(tem.y, weaBus.TDryBul) annotation (Line(points={{-29,60},{
          -12,60},{-12,56},{68,56},{68,-2}},
                        color={0,0,127}), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  connect(pressure.y, weaBus.pAtm) annotation (Line(points={{-27,-28},
          {12,-28},{12,-32},{68,-32},{68,-2}},
                              color={0,0,127}), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));

  connect(relativeHum.y, weaBus.relHum) annotation (Line(points={{-23,-56},
          {18,-56},{18,-58},{68,-58},{68,-2}},
                              color={0,0,127}), Text(
      string="%second",
      index=1,
      extent={{6,3},{6,3}}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={Bitmap(extent={{-56,62},{64,
              -74}}, fileName="")}));
end myWea;
