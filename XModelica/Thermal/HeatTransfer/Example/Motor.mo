within Exergy.XModelica.Thermal.HeatTransfer.Example;
model Motor "Second order thermal model of a motor"
  extends Modelica.Icons.Example;
  parameter Modelica.SIunits.Temperature TAmb(displayUnit="degC") = 293.15
    "Ambient temperature";

  Modelica.Blocks.Sources.CombiTimeTable lossTable(extrapolation=Modelica.
        Blocks.Types.Extrapolation.Periodic, table=[0,100,500; 360,100,500;
        360,1000,500; 600,1000,500])
                            annotation (Placement(transformation(
        origin={-40,70},
        extent={{-10,-10},{10,10}},
        rotation=270)));
   Sources.PrescribedHeatFlow windingLosses(T_ref=368.15, alpha=
        3.03E-3) annotation (Placement(transformation(
        origin={-80,10},
        extent={{-10,-10},{10,10}},
        rotation=270)));
   Components.HeatCapacitor winding(C=2500, T(start=TAmb, fixed=true))
    annotation (Placement(transformation(extent={{-90,-20},{-70,-40}},
          rotation=0)));
   Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor Twinding annotation (
      Placement(transformation(
        origin={-60,-50},
        extent={{-10,-10},{10,10}},
        rotation=270)));
   Components.ThermalConductor winding2core(G=10) annotation (
      Placement(transformation(extent={{-50,-20},{-30,0}}, rotation=0)));
   Sources.PrescribedHeatFlow coreLosses annotation (Placement(
        transformation(
        origin={0,10},
        extent={{-10,-10},{10,10}},
        rotation=270)));
   Components.HeatCapacitor core(C=25000, T(start=TAmb, fixed=true))
    annotation (Placement(transformation(extent={{-10,-20},{10,-40}},
          rotation=0)));
  Modelica.Thermal.HeatTransfer.Celsius.TemperatureSensor Tcore annotation (
      Placement(transformation(
        origin={-20,-50},
        extent={{-10,-10},{10,10}},
        rotation=270)));
  Modelica.Blocks.Sources.Constant convectionConstant(k=25)
    annotation (Placement(transformation(
        origin={40,30},
        extent={{-10,-10},{10,10}},
        rotation=270)));
   Components.Convection convection annotation (Placement(
        transformation(extent={{30,-20},{50,0}}, rotation=0)));
   Sources.FixedTemperature environment(T=TAmb) annotation (Placement(
        transformation(
        origin={80,-10},
        extent={{-10,-10},{10,10}},
        rotation=180)));
      //    Exergy.Utilities.ViewObject viewObject(nVolume=1);
  inner XExergy.Utilities.System system annotation (Placement(
        transformation(extent={{-100,72},{-80,92}})));
  XExergy.Utilities.ViewRoute viewRoute(nSubSystem=8)
    annotation (Placement(transformation(extent={{64,80},{84,100}})));
  XExergy.Utilities.ViewPort viewPort annotation (Placement(
        transformation(extent={{90,80},{110,100}})));
  XExergy.Utilities.ViewObject viewObject1(nVolume=1)
    annotation (Placement(transformation(extent={{64,58},{84,78}})));
equation
  connect(windingLosses.port, winding.port)  annotation (Line(points={{-80,0},
          {-80,-20}},    color={191,0,0}));
  connect(coreLosses.port, core.port)  annotation (Line(points={{0,0},{0,
          -10},{0,-20}},                                       color={191,0,
          0}));
  connect(winding.port, winding2core.port_a)
                                   annotation (Line(points={{-80,-20},{-80,
          -10},{-50,-10}}, color={191,0,0}));
  connect(winding2core.port_b, core.port)
                                annotation (Line(points={{-30,-10},{0,-10},
          {0,-20}}, color={191,0,0}));
  connect(winding.port, Twinding.port)  annotation (Line(points={{-80,-20},
          {-80,-10},{-60,-10},{-60,-40}}, color={191,0,0}));
  connect(core.port, Tcore.port)  annotation (Line(points={{0,-20},{0,-10},
          {-20,-10},{-20,-40}}, color={191,0,0}));
  connect(winding2core.port_b, convection.solid)
                                      annotation (Line(points={{-30,-10},{
          30,-10}}, color={191,0,0}));
  connect(convection.fluid, environment.port) annotation (Line(points={{50,-10},
          {60,-10},{70,-10}},               color={191,0,0}));
  connect(convectionConstant.y, convection.Gc)
    annotation (Line(points={{40,19},{40,0}}, color={0,0,127}));
  connect(lossTable.y[1], windingLosses.Q_flow) annotation (Line(points={{-40,59},
          {-40,40},{-80,40},{-80,20}},         color={0,0,127}));
  connect(lossTable.y[2], coreLosses.Q_flow) annotation (Line(points={{-40,59},
          {-40,40},{0,40},{0,20}},                             color={0,0,
          127}));
          /*
  viewObject.volume[1].M = sum(winding2core.viewObject.volume.M) +
    sum(winding.viewObject.volume.M) + sum(core.viewObject.volume.M)
    +sum(environment.viewObject.volume.M)+sum(convection.viewObject.volume.M)
    +sum(windingLosses.viewObject.volume.M)+sum(coreLosses.viewObject.volume.M);
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
  connect(windingLosses.viewOutput, viewRoute.viewOutput[1]) annotation (Line(
      points={{-70.2,0.4},{-70,0.4},{-70,87.175},{63.9,87.175}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(coreLosses.viewOutput, viewRoute.viewOutput[2]) annotation (Line(
      points={{9.8,0.4},{6,0.4},{6,87.925},{63.9,87.925}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(winding2core.viewOutput, viewRoute.viewOutput[3]) annotation (Line(
      points={{-30.4,-0.2},{-30.4,18},{-30.4,88.675},{63.9,88.675}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(core.viewOutput, viewRoute.viewOutput[4]) annotation (Line(
      points={{9.6,-39.8},{9.6,-32},{9.6,89.425},{63.9,89.425}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(winding.viewOutput, viewRoute.viewOutput[5]) annotation (Line(
      points={{-70.4,-39.8},{-64,-39.8},{-64,84},{63.9,84},{63.9,
          90.175}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(convection.viewOutput, viewRoute.viewOutput[6]) annotation (Line(
      points={{49.6,-0.2},{63.9,-0.2},{63.9,90.925}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(environment.viewOutput, viewRoute.viewOutput[7]) annotation (Line(
      points={{70.4,-19.8},{70.4,-2},{70.4,64},{42,64},{42,91.675},{
          63.9,91.675}},
      color={28,108,200},
      pattern=LinePattern.Dot));

  connect(viewRoute.viewTotal, viewPort) annotation (Line(
      points={{83.8,89.9},{89.9,89.9},{89.9,90},{100,90}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  connect(viewObject1.viewOutput, viewRoute.viewTotal) annotation (Line(
      points={{84,68},{84,89.9},{83.8,89.9}},
      color={28,108,200},
      pattern=LinePattern.Dot));
  annotation (Documentation(info="<HTML>
<p>
This example contains a simple second order thermal model of a motor.
The periodic power losses are described by table \"lossTable\":
</p>
<table>
<tr><td valign=\"top\">time</td><td valign=\"top\">winding losses</td><td valign=\"top\">core losses</td></tr>
<tr><td valign=\"top\">   0</td><td valign=\"top\">           100</td><td valign=\"top\">        500</td></tr>
<tr><td valign=\"top\"> 360</td><td valign=\"top\">           100</td><td valign=\"top\">        500</td></tr>
<tr><td valign=\"top\"> 360</td><td valign=\"top\">          1000</td><td valign=\"top\">        500</td></tr>
<tr><td valign=\"top\"> 600</td><td valign=\"top\">          1000</td><td valign=\"top\">        500</td></tr>
</table>
<p>
Since constant speed is assumed, the core losses keep constant
whereas the winding losses are low for 6 minutes (no-load) and high for 4 minutes (over load).
</p>
<p>
The winding losses are corrected by (1 + alpha*(T - T_ref)) because the winding's resistance is temperature dependent whereas the core losses are kept constant (alpha = 0).
</p>
<p>
The power dissipation to the environment is approximated by heat flow through
a thermal conductance between winding and core,
partially storage of the heat in the winding's heat capacity
as well as the core's heat capacity and finally by forced convection to the environment.<br>
Since constant speed is assumed, the convective conductance keeps constant.<br>
Using Modelica.Thermal.FluidHeatFlow it would be possible to model the coolant air flow, too
(instead of simple dissipation to a constant ambient's temperature).
</p>
<p>
Simulate for 7200 s; plot Twinding.T and Tcore.T.
</p>
</HTML>"),
    experiment(StopTime=720, Interval=0.01),
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}})));
end Motor;
