within Exergy.XThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor;
model Test_nozzle
  import ThermoCycle;

  ThermoCycle.Components.FluidFlow.Reservoirs.SourceMdot sourceMdot(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
    h_0=480000,
    Mdot_0=0.005,
    p=1800000,
    UseT=true,
    T_0=423.15)
    annotation (Placement(transformation(extent={{-140,-24},{-104,12}})));
  ThermoCycle.Components.FluidFlow.Reservoirs.SinkP sinkP(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
    h=430000,
    p0=600000)
    annotation (Placement(transformation(extent={{0,-22},{20,-2}})));
  Nozzle         fuite_interne2(
  redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
  Use_gamma=false,
    P_su_start=1800000,
    T_su_start=423.15,
    P_ex_start=600000)
    annotation (Placement(transformation(extent={{-88,-28},{-46,18}})));
  Modelica.Blocks.Sources.Ramp ramp(
    duration=5,
    startTime=0,
    height=0.01,
    offset=0.0001)
    annotation (Placement(transformation(extent={{-174,-14},{-154,6}})));
equation
  connect(ramp.y, sourceMdot.in_Mdot) annotation (Line(
      points={{-153,-4},{-144,-4},{-144,28},{-132.8,28},{-132.8,4.8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sourceMdot.flangeB, fuite_interne2.su) annotation (Line(
      points={{-105.8,-6},{-98,-6},{-98,-4.54},{-88.42,-4.54}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(fuite_interne2.ex, sinkP.flangeB) annotation (Line(
      points={{-46,-4.54},{-22,-4.54},{-22,-12},{1.6,-12}},
      color={0,0,255},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),      graphics),
    experiment(__Dymola_NumberOfIntervals=10),
    __Dymola_experimentSetupOutput);
end Test_nozzle;
