within Exergy.XClaRa.Components.Utilities.Blocks.Check;
model test_1_LimPID
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
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb50;
  LimPID         PID(
    k=10,
    Tau_in=1e-3,
    Tau_out=1e-3,
    y_start=4,
    xi_start=3,
    xd_start=2,
    initType=Modelica.Blocks.Types.InitPID.InitialOutput,
    controllerType=Modelica.Blocks.Types.SimpleController.PI,
    y_max=100,
    Tau_i=0.1,
    y_inactive=0)
           annotation (Placement(transformation(extent={{-38,-20},{-18,0}})));
  Modelica.Blocks.Sources.RealExpression realExpression(y=1)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
  Modelica.Blocks.Sources.Sine sine(
    amplitude=0.2,
    freqHz=0.1,
    offset=0) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={72,-62})));
  Modelica.Blocks.Continuous.FirstOrder firstOrder(
    initType=Modelica.Blocks.Types.Init.InitialOutput,
    T=1,
    y_start=0.3)
         annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={-6,-56})));
  Modelica.Blocks.Math.Add add annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=180,
        origin={22,-56})));
  Modelica.Blocks.Sources.BooleanExpression booleanExpression(y=time > 10)
    annotation (Placement(transformation(extent={{-80,-8},{-60,-28}})));
equation
  connect(realExpression.y, PID.u_s) annotation (Line(
      points={{-59,-10},{-40,-10}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(add.y, firstOrder.u) annotation (Line(
      points={{11,-56},{6,-56}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sine.y, add.u1) annotation (Line(
      points={{61,-62},{34,-62}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(PID.y, add.u2) annotation (Line(
      points={{-17,-10},{52,-10},{52,-50},{34,-50}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(firstOrder.y, PID.u_m) annotation (Line(
      points={{-17,-56},{-27.9,-56},{-27.9,-22}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (
    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{
            100,100}}),
            graphics={            Text(
          extent={{-94,98},{104,58}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
PURPOSE:

______________________________________________________________________________________________
"),                    Text(
          extent={{-134,102},{66,82}},
          lineColor={0,128,0},
          fontSize=31,
          textString="TESTED -- YYYY-MM-DD //XX"),Text(
          extent={{-94,58},{70,44}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          textString="______________________________________________________________________________________________________________
Remarks: 
______________________________________________________________________________________________________________
",        fontSize=8),Text(
          extent={{-94,72},{106,54}},
          lineColor={0,128,0},
          horizontalAlignment=TextAlignment.Left,
          fontSize=10,
          textString="______________________________________________________________________________________________
Scenario:  

______________________________________________________________________________________________
")}),
    experiment(StopTime=50),
    __Dymola_experimentSetupOutput(equdistant=false));
end test_1_LimPID;
