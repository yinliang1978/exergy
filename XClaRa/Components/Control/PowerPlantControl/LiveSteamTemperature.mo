within Exergy.XClaRa.Components.Control.PowerPlantControl;
model LiveSteamTemperature
  "A simple controller for the live steam temperature based on Strauss: \"Kraftwerkstechnik\", 5th edition, 2006."
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

  Utilities.Blocks.LimPID PID2(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    sign=1,
    u_ref=T_a2_ref,
    y_ref=T_e2_ref,
    k=k_PID2,
    Tau_i=Tau_i_PID2,
    y_max=1000,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    y_start=536,
    perUnitConversion=false) annotation (Placement(transformation(
        extent={{10,9},{-10,-9}},
        rotation=90,
        origin={7,68})));
  Utilities.Blocks.LimPID feedback2(
    controllerType=Modelica.Blocks.Types.SimpleController.P,
    perUnitConversion=false,
    y_min=0,
    y_max=1,
    k=k_P2,
    sign=-1,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    y_start=0.9) annotation (Placement(transformation(
        extent={{10,9.5},{-10,-9.5}},
        rotation=90,
        origin={7.5,36})));
  Modelica.Blocks.Math.Add add2 annotation (Placement(transformation(
        extent={{10,-10},{-10,10}},
        rotation=90,
        origin={56,36})));
  Modelica.Blocks.Interfaces.RealOutput opening2 "valve opening of injector 2"
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-110,20})));
public
  ClaRa.Basics.Interfaces.Bus MeasurementValues annotation (Placement(
        transformation(extent={{-20,78},{20,118}}),
        iconTransformation(extent={{-8,88},{12,110}})));
public
  ClaRa.Basics.Interfaces.Bus SetValues annotation (Placement(
        transformation(extent={{30,78},{70,118}}), iconTransformation(
          extent={{40,90},{60,110}})));
  Utilities.Blocks.LimPID PID1(
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    u_ref=Delta_T_2_ref,
    y_ref=T_e1_ref,
    k=k_PID1,
    Tau_i=Tau_i_PID1,
    y_max=1000,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    y_start=492,
    sign=+1,
    perUnitConversion=false) annotation (Placement(transformation(
        extent={{10,9},{-10,-9}},
        rotation=90,
        origin={-2,-28})));
  Utilities.Blocks.LimPID P1(
    controllerType=Modelica.Blocks.Types.SimpleController.P,
    sign=-1,
    perUnitConversion=false,
    y_max=1,
    y_min=0,
    k=k_P1,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    y_start=0.5) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={-2,-68})));
  Modelica.Blocks.Interfaces.RealOutput opening1 "opening of  injector 1"
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-110,-80})));
  parameter Real T_a2_ref=873.15 "Reference value for controlled variable";
  parameter Real k_PID2=1 "Gain of PID2";
  parameter Modelica.SIunits.Time Tau_i_PID2=0.5
    "Integration time constant of PID2";
  parameter Real T_e2_ref=1 "Reference value injector 2 outlet temperature";
  parameter Real k_P2=1 "Gain value multiplied with input signal";
  parameter Real Delta_T_2_ref=20
    "Reference value for Temperature difference over injector 2";
  parameter Real T_e1_ref=1 "Reference value for injector 1 outlet temperature";
  parameter Real k_PID1=1 "Gain of controller";
  parameter Modelica.SIunits.Time Tau_i_PID1=0.5
    "Time constant of Integrator block";
  parameter Real k_P1=1 "Gain of controller";
equation
  connect(PID2.y, add2.u1) annotation (Line(
      points={{7,57},{50,57},{50,48}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(MeasurementValues.T_a2, PID2.u_m) annotation (Line(
      points={{0,98},{-72,98},{-72,68},{-3.8,68}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(SetValues.T_a2_set, PID2.u_s) annotation (Line(
      points={{50,98},{44.5,98},{44.5,80},{7,80}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(SetValues.Delta_T2_set, add2.u2) annotation (Line(
      points={{50,98},{62,98},{62,48}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(MeasurementValues.T_a1, PID1.u_m) annotation (Line(
      points={{0,98},{-72,98},{-72,-28},{-12.8,-28}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(PID1.u_s, add2.y) annotation (Line(
      points={{-2,-16},{56,-16},{56,25}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(P1.y, opening1) annotation (Line(
      points={{-2,-79},{-2,-80},{-110,-80}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PID1.y, P1.u_s) annotation (Line(
      points={{-2,-39},{-2,-56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(feedback2.y, opening2) annotation (Line(
      points={{7.5,25},{10,25},{10,20},{-110,20}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PID2.y, feedback2.u_s) annotation (Line(
      points={{7,57},{7,52.5},{7.5,52.5},{7.5,48}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(MeasurementValues.T_e2, feedback2.u_m) annotation (Line(
      points={{0,98},{-72,98},{-72,36},{-3.9,36}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  connect(MeasurementValues.T_e1, P1.u_m) annotation (Line(
      points={{0,98},{-72,98},{-72,-68},{-14,-68}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None), Text(
      string="%first",
      index=-1,
      extent={{-6,3},{-6,3}}));
  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-100,100},{0,2},{-100,-100}},
          color={0,0,255},
          smooth=Smooth.None),
        Text(
          extent={{2,42},{84,-40}},
          lineColor={0,0,255},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          textString="T")}));
end LiveSteamTemperature;
