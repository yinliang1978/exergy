within Exergy.XThermoCycle.Components.Units.Tanks.HeatStorageWaterHeater;
model Test_heat_storage

  Exergy.XThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot sourceMdot(
    h_0=2.8E5,
    redeclare package Medium = ThermoCycle.Media.StandardWater,
    UseT=true,
    Mdot_0=0.2,
    p=100000,
    T_0=283.15) annotation (Placement(transformation(extent={{-82,-88},
            {-50,-56}})));
  Exergy.XThermoCycle.Components.FluidFlow.Reservoirs.SinkP sinkP(h=2E5,
      redeclare package Medium = ThermoCycle.Media.StandardWater)
    annotation (Placement(transformation(extent={{28,22},{52,44}})));
  Exergy.XThermoCycle.Components.HeatFlow.Sources.Source_T_cell source_T
    annotation (Placement(transformation(extent={{-60,-12},{-40,8}})));
  Modelica.Blocks.Sources.Constant const(k=273.15 + 25)
    annotation (Placement(transformation(extent={{-76,24},{-56,44}})));

  Exergy.XThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot sourceMdot1(
    h_0=2.8E5,
    UseT=true,
    redeclare package Medium = ThermoCycle.Media.StandardWater,
    Mdot_0=1,
    p=100000,
    T_0=363.15) annotation (Placement(transformation(extent={{82,-28},
            {50,4}})));
  Exergy.XThermoCycle.Components.FluidFlow.Reservoirs.SinkP sinkP1(h=2E5,
      redeclare package Medium = ThermoCycle.Media.StandardWater)
    annotation (Placement(transformation(extent={{38,-72},{62,-50}})));

  Heat_storage_hx heat_storage_hx(
    N=25,
    V_tank=0.15,
    h1=0.2,
    h2=0.8,
    Unom_hx=2000)
    annotation (Placement(transformation(extent={{-36,-66},{36,10}})));
  Modelica.Blocks.Sources.Step step(
    offset=273.15 + 90,
    height=-40,
    startTime=100)
    annotation (Placement(transformation(extent={{26,66},{46,86}})));
equation

  connect(const.y, source_T.Temperature) annotation (Line(
      points={{-55,34},{-44,34},{-44,3},{-49,3}},
      color={0,0,127},
      smooth=Smooth.None));

  connect(sourceMdot1.flangeB, heat_storage_hx.SecondaryFluid_su) annotation (
      Line(
      points={{51.6,-12},{32,-12},{32,-15.08},{14.4,-15.08}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(heat_storage_hx.SecondaryFluid_ex, sinkP1.flangeB) annotation (Line(
      points={{14.4,-39.4},{27.2,-39.4},{27.2,-61},{39.92,-61}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(sinkP.flangeB, heat_storage_hx.MainFluid_ex) annotation (Line(
      points={{29.92,33},{0,33},{0,4.68}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(sourceMdot.flangeB, heat_storage_hx.MainFluid_su) annotation (Line(
      points={{-51.6,-72},{-32,-72},{-32,-58.02},{-13.32,-58.02}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(source_T.ThermalPortCell, heat_storage_hx.Wall_ext) annotation (
      Line(
      points={{-49.1,-5.1},{-49.1,-26.86},{-13.32,-26.86}},
      color={255,0,0},
      smooth=Smooth.None));
  connect(step.y, sourceMdot1.in_T) annotation (Line(
      points={{47,76},{66.32,76},{66.32,-2.4}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Line(
      points={{-0.4,-12.26},{-0.4,-9.445},{-0.08,-9.445},{-0.08,-6.63},{2.32,
          -6.63},{2.32,-1.5}},
      color={255,0,0},
      smooth=Smooth.None), Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics),
    experiment(StopTime=1000),
    __Dymola_experimentSetupOutput);
end Test_heat_storage;
