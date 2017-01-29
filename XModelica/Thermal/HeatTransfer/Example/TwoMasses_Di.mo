within Exergy.XModelica.Thermal.HeatTransfer.Example;
model TwoMasses_Di "Simple conductor demo"
  extends Modelica.Icons.Example;

  Components.HeatCapacitor_Di mass1(C=15, T(fixed=true, start=373.15))
    annotation (Placement(transformation(extent={{-156,6},{-88,74}})));
  Components.HeatCapacitor_Di mass2(C=15, T(fixed=true, start=273.15))
    annotation (Placement(transformation(extent={{2,16},{64,78}})));
  Components.ThermalConductor conductor(G=10)
    annotation (Placement(transformation(extent={{-80,-46},{-10,24}})));
  Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor Tsensor1
    annotation (Placement(transformation(extent={{-108,-86},{-70,-48}})));
  Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor Tsensor2
    annotation (Placement(transformation(extent={{22,-86},{-18,-46}})));

  // attention here vVolume must equate to 1
  Utilities.ViewObject viewObject(nEnergy={0,1,0,0});
  inner Exergy.Utilities.RefEnv   refEnv;
  Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{106,88},{126,108}})));
  Utilities.ViewRoute viewRoute(nSubSystem=3)
    annotation (Placement(transformation(extent={{76,86},{96,106}})));
equation
  connect(mass1.port, conductor.port_a) annotation (Line(points={{-122,6},{-122,
          6},{-122,-11},{-80,-11}}, color={191,0,0}));
  connect(Tsensor1.port, conductor.port_a) annotation (Line(points={{-108,-67},{
          -122,-67},{-122,-11},{-80,-11}}, color={191,0,0}));
/*
  viewObject.volume[1].M = sum(mass1.viewObject.volume.M) +
    sum(mass2.viewObject.volume.M) + sum(conductor.viewObject.volume.M);
  viewObject.volume[1].E = sum(mass1.viewObject.volume.E) +
    sum(mass2.viewObject.volume.E) + sum(conductor.viewObject.volume.E);
  viewObject.volume[1].S = sum(mass1.viewObject.volume.S) +
    sum(mass2.viewObject.volume.S) + sum(conductor.viewObject.volume.S);
  viewObject.volume[1].V = 0;

  viewObject.volume[1].E_0 = sum(mass1.viewObject.volume.E_0)
     + sum(mass2.viewObject.volume.E_0) + sum(
    conductor.viewObject.volume.E_0);
  viewObject.volume[1].S_0 = sum(mass1.viewObject.volume.S_0)
     + sum(mass2.viewObject.volume.S_0) + sum(
    conductor.viewObject.volume.S_0);
    viewObject.volume[1].V_0 = 0;
    
    */
    connect(viewObject.viewOutput,viewRoute.viewTotal);

  connect(Tsensor2.port, mass2.port) annotation (Line(points={{22,-66},{32,-66},
          {32,-8},{32,16},{33,16}}, color={191,0,0}));
  connect(conductor.port_b, mass2.port) annotation (Line(points={{-10,-11},{32,-11},
          {32,16},{33,16}}, color={191,0,0}));
  connect(viewRoute.viewOutput[1], mass1.viewOutput) annotation (Line(
      points={{75.9,93.8},{66,93.8},{-89.36,93.8},{-89.36,73.32}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(mass2.viewOutput, viewRoute.viewOutput[2]) annotation (Line(
      points={{62.76,77.38},{62.76,95.8},{75.9,95.8}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(conductor.viewOutput, viewRoute.viewOutput[3]) annotation (Line(
      points={{-11.4,23.3},{-11.4,97.8},{75.9,97.8}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(viewOutput, viewRoute.viewTotal) annotation (Line(
      points={{116,98},{95.8,98},{95.8,95.9}},
      color={28,108,200},
      pattern=LinePattern.Dot));
    annotation (Placement(transformation(extent={{82,40},{102,60}})),
              Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-160,-100},{120,100}})), Icon(coordinateSystem(extent={{-160,
            -100},{120,100}})));
end TwoMasses_Di;
