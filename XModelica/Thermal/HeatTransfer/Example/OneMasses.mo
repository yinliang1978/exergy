within Exergy.XModelica.Thermal.HeatTransfer.Example;
model OneMasses "Simple conductor demo"
  extends Modelica.Icons.Example;

  Components.HeatCapacitor mass2(C=15, T(fixed=true, start=273.15))
    annotation (Placement(transformation(extent={{2,16},{64,78}})));
  Components.ThermalConductor conductor(G=10)
    annotation (Placement(transformation(extent={{-80,-46},{-10,24}})));
  Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor Tsensor1
    annotation (Placement(transformation(extent={{-108,-86},{-70,-48}})));
  Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor Tsensor2
    annotation (Placement(transformation(extent={{22,-86},{-18,-46}})));

//  Exergy.Utilities.ViewObject viewObject(nVolume=1);
  inner XExergy.Utilities.System system annotation (Placement(
        transformation(extent={{-158,76},{-138,96}})));
  Sources.FixedTemperature mass1(T=373.15)
    annotation (Placement(transformation(extent={{-23,-23},{23,23}},
        rotation=-90,
        origin={-121,37})));
  XExergy.Utilities.ViewRoute viewRoute(nSubSystem=3) annotation (
      Placement(transformation(extent={{86,82},{106,102}})));
  XExergy.Utilities.ViewPort viewPort annotation (Placement(
        transformation(extent={{110,82},{130,102}})));
  XExergy.Utilities.ViewObject viewObject(nVolume=1)
    annotation (Placement(transformation(extent={{84,54},{104,74}})));
equation
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
      //heat
 // viewObject.heat[1].E_flow = conductor.port_a.Q_flow;
 // viewObject.heat[1].T = conductor.port_a.T;
/*
  //volume
      viewObject.volume[1].M =
    sum(mass2.viewObject.volume.M) + sum(conductor.viewObject.volume.M);
  viewObject.volume[1].E =
    sum(mass2.viewObject.volume.E) + sum(conductor.viewObject.volume.E);
  viewObject.volume[1].S =
    sum(mass2.viewObject.volume.S) + sum(conductor.viewObject.volume.S);
  viewObject.volume[1].V = 0;

  viewObject.volume[1].E_0 =  sum(mass2.viewObject.volume.E_0) + sum(
    conductor.viewObject.volume.E_0);
  viewObject.volume[1].S_0 = sum(mass2.viewObject.volume.S_0) + sum(
    conductor.viewObject.volume.S_0);
    viewObject.volume[1].V_0 = 0;
*/
  connect(Tsensor2.port, mass2.port) annotation (Line(points={{22,-66},{32,-66},
          {32,-8},{32,16},{33,16}}, color={191,0,0}));
  connect(conductor.port_b, mass2.port) annotation (Line(points={{-10,-11},{32,-11},
          {32,16},{33,16}}, color={191,0,0}));
  connect(mass1.port, conductor.port_a) annotation (Line(points={{-121,14},{-121,
          14},{-122,14},{-122,-11},{-80,-11}},                     color={191,0,
          0}));
  connect(viewRoute.viewTotal, viewPort) annotation (Line(points={{105.8,91.9},{
          112.9,91.9},{112.9,92},{120,92}}, color={28,108,200},
      pattern=LinePattern.Dot));
  connect(mass1.viewOutput, viewRoute.viewOutput[1]) annotation (Line(points={{-98.46,
          14.92},{-94,14.92},{-94,89.8},{85.9,89.8}}, color={28,108,
          200},
      pattern=LinePattern.Dot));
  connect(mass2.viewOutput, viewRoute.viewOutput[2]) annotation (Line(points={{62.76,
          77.38},{62.76,94},{85.9,94},{85.9,91.8}}, color={28,108,200},
      pattern=LinePattern.Dot));

  connect(conductor.viewOutput, viewRoute.viewOutput[3]) annotation (Line(
        points={{-11.4,23.3},{-11.4,93.8},{85.9,93.8}}, color={28,108,
          200},
      pattern=LinePattern.Dot));
  connect(viewObject.viewOutput, viewRoute.viewTotal) annotation (
      Line(
      points={{104,64},{104,91.9},{105.8,91.9}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-160,-100},{120,100}})), Icon(coordinateSystem(extent={{-160,
            -100},{120,100}})));
end OneMasses;
