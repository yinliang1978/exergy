within Exergy.XModelica.Thermal.HeatTransfer.Example;
model TwoMassesTest "Simple conductor demo"
  extends Modelica.Icons.Example;

  Components.HeatCapacitor mass1(C=15, T(fixed=true, start=373.15))
    annotation (Placement(transformation(extent={{-156,6},{-88,74}})));
  Components.HeatCapacitor mass2(C=15, T(fixed=true, start=273.15))
    annotation (Placement(transformation(extent={{2,16},{64,78}})));
  Components.ThermalConductor conductor(G=10)
    annotation (Placement(transformation(extent={{-80,-46},{-10,24}})));

  // attention here vVolume must equate to 1

  inner Exergy.Utilities.RefEnv   refEnv;

  Utilities.ViewRoute viewRoute(nSubSystem=3)
    annotation (Placement(transformation(extent={{88,84},{108,104}})));
  Utilities.ViewPort viewPort
    annotation (Placement(transformation(extent={{112,84},{132,104}})));
  Utilities.ViewObject viewObject
    annotation (Placement(transformation(extent={{70,50},{90,70}})));
equation
  connect(mass1.port, conductor.port_a) annotation (Line(points={{-122,6},{-122,
          6},{-122,-11},{-80,-11}}, color={191,0,0}));
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

  connect(conductor.port_b, mass2.port) annotation (Line(points={{-10,-11},{32,-11},
          {32,16},{33,16}}, color={191,0,0}));
  connect(mass1.viewOutput, viewRoute.viewOutput[1]) annotation (Line(points={{-89.36,
          73.32},{-89.36,91.8},{87.9,91.8}}, color={28,108,200}));
  connect(mass2.viewOutput, viewRoute.viewOutput[2]) annotation (Line(points={{62.76,
          77.38},{62.76,93.8},{87.9,93.8}}, color={28,108,200}));
  connect(viewRoute.viewTotal, viewPort) annotation (Line(points={{107.8,93.9},{
          116.9,93.9},{116.9,94},{122,94}}, color={28,108,200}));
  connect(conductor.viewOutput, viewRoute.viewOutput[3]) annotation (Line(
        points={{-11.4,23.3},{-11.4,95.8},{87.9,95.8}}, color={28,108,200}));
  connect(viewObject.viewOutput, viewPort)
    annotation (Line(points={{90,60},{122,60},{122,94}}, color={28,108,200}));
    annotation (Placement(transformation(extent={{82,40},{102,60}})),
              Diagram(coordinateSystem(preserveAspectRatio=false,
          extent={{-160,-100},{120,100}})), Icon(coordinateSystem(extent={{-160,
            -100},{120,100}})));
end TwoMassesTest;
