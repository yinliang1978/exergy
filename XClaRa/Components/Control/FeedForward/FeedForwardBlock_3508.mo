within Exergy.XClaRa.Components.Control.FeedForward;
model FeedForwardBlock_3508
  "feed forward for coal mass flow and turbine valve aperture"
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                        //
//                                                                           //
// Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
// Copyright © 2013-2016, DYNCAP/DYNSTART research team.                     //
//___________________________________________________________________________//
// DYNCAP and DYNSTART are research projects supported by the German Federal //
// Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
// The research team consists of the following project partners:             //
// Institute of Energy Systems (Hamburg University of Technology),           //
// Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
// TLK-Thermo GmbH (Braunschweig, Germany),                                  //
// XRG Simulation GmbH (Hamburg, Germany).                                   //
//___________________________________________________________________________//

  parameter Real CL_Valve_[:,:]=[0,0; 1, 1]
    "Characteristics of the turbine valve" annotation(Dialog(group="Part Load Definition"));

  parameter Real k_FuelOvershoot=100 "Gain of fuel overshoot";
  parameter Real eta_gen= 0.98 "Efficiency of generator";

  Modelica.Blocks.Interfaces.RealInput derP_max_
    "Maximum overall power gradient" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=180,
        origin={120,-20})));
  Modelica.Blocks.Interfaces.RealInput derP_StG_
    "Maximum  power gradient due to steam generator restrictions" annotation (
      Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=180,
        origin={120,-50})));
  Modelica.Blocks.Interfaces.RealInput derP_T_
    "Maximum  power gradient due to turbine restrictions" annotation (Placement(
        transformation(
        extent={{-20,-20},{20,20}},
        rotation=180,
        origin={120,-80})));
  Modelica.Blocks.Interfaces.RealInput P_G_target_
    "Target value generator power" annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=270,
        origin={-40,120})));
protected
  Modelica.Blocks.Math.Min min annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={76,32})));
public
  Modelica.Blocks.Math.Min min1 annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={40,26})));
  Modelica.Blocks.Nonlinear.VariableLimiter variableLimiter annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-76,68})));
  Utilities.Blocks.VariableGradientLimiter variableGradientLimiter annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-76,10})));
  Modelica.Blocks.Math.Gain gain(k=-1) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-42,40})));
  Modelica.Blocks.Interfaces.RealInput P_max_ "Maximum generator power"
                                     annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=180,
        origin={120,80})));
  Modelica.Blocks.Interfaces.RealInput P_min_ "Minimum generator power"
                                     annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=180,
        origin={120,50})));
  Modelica.Blocks.Math.Gain FuelFeedForward(k=1/eta_gen)
                                                 annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-76,-36})));
  Modelica.Blocks.Continuous.Derivative FuelOvershoot(initType=Modelica.Blocks.Types.Init.InitialOutput,
      k=k_FuelOvershoot)                           annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-38,-32})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-70,-78})));
  Modelica.Blocks.Interfaces.RealOutput QF_FF_ "FiringPower feed forward"
                                  annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-40,-110})));
  Modelica.Blocks.Tables.CombiTable1D turbineValveOpeneing(table=CL_Valve_,
      columns={2}) "load dependend turbine valve opening"
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={16,-38})));
  Modelica.Blocks.Interfaces.RealOutput y_Turbine_
    "Feed forward value for turbine valve openeing" annotation (Placement(
        transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={20,-110})));
  Modelica.Blocks.Interfaces.RealOutput P_G_set_
    "Connector of Real output signal" annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-80,-110})));
equation
  connect(derP_max_, min.u2) annotation (Line(
      points={{120,-20},{102,-20},{102,38},{88,38}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(derP_StG_, min.u1) annotation (Line(
      points={{120,-50},{98,-50},{98,26},{88,26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(min.y, min1.u2) annotation (Line(
      points={{65,32},{65,31},{52,31},{52,32}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(derP_T_, min1.u1) annotation (Line(
      points={{120,-80},{94,-80},{94,20},{52,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P_G_target_, variableLimiter.u) annotation (Line(
      points={{-40,120},{-40,80},{-76,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(min1.y, variableGradientLimiter.maxGrad) annotation (Line(
      points={{29,26},{-68,26},{-68,22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableLimiter.y, variableGradientLimiter.u) annotation (Line(
      points={{-76,57},{-76,22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.u, min1.y) annotation (Line(
      points={{-30,40},{22.5,40},{22.5,26},{29,26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(gain.y, variableGradientLimiter.minGrad) annotation (Line(
      points={{-53,40},{-84.5,40},{-84.5,22},{-84,22}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableLimiter.limit1, P_max_) annotation (Line(
      points={{-68,80},{120,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P_min_, variableLimiter.limit2) annotation (Line(
      points={{120,50},{68,50},{68,90},{-84,90},{-84,80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableGradientLimiter.y, FuelFeedForward.u) annotation (Line(
      points={{-76,-1},{-76,-24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(FuelOvershoot.u, variableGradientLimiter.y)
                                                   annotation (Line(
      points={{-38,-20},{-38,-10},{-76,-10},{-76,-1}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u2, FuelFeedForward.y) annotation (Line(
      points={{-76,-66},{-76,-47}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.u1, FuelOvershoot.y)
                                annotation (Line(
      points={{-64,-66},{-64,-62},{-38,-62},{-38,-43}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, QF_FF_)     annotation (Line(
      points={{-70,-89},{-40,-89},{-40,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableGradientLimiter.y, turbineValveOpeneing.u[1]) annotation (
      Line(
      points={{-76,-1},{-76,-10},{16,-10},{16,-26}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(turbineValveOpeneing.y[1], y_Turbine_) annotation (Line(
      points={{16,-49},{20,-49},{20,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(variableGradientLimiter.y, P_G_set_) annotation (Line(
      points={{-76,-1},{-96,-1},{-96,-94},{-80,-94},{-80,-110}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-94,-40},{-62,-40},{-36,-40},{2,-40},{-6,-16},{-40,6},{-64,-18},
              {-78,-10},{-94,-40}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-100,58},{102,34}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="FeedForwardBlock")}));
end FeedForwardBlock_3508;
