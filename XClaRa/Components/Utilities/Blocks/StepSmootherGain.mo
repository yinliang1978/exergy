within Exergy.XClaRa.Components.Utilities.Blocks;
block StepSmootherGain "Smoothly activate and deactivate a Real signal"
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                            //
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
  extends Modelica.Blocks.Interfaces.SISO;

  Modelica.Blocks.Interfaces.RealInput x
    "y = u for x > func and y=0 for x < noFunc"                                      annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-120}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={0,-120})));
  Modelica.Blocks.Interfaces.RealInput noFunc "y=0 for x < noFunc"  annotation (Placement(transformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={60,-120}), iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-70,-120})));
  Modelica.Blocks.Interfaces.RealInput func "y = u for x > func"  annotation (Placement(transformation(extent={{-20,-20},{20,20}},
        rotation=90,
        origin={-60,-120}),                                               iconTransformation(
        extent={{-20,-20},{20,20}},
        rotation=90,
        origin={68,-120})));

equation
  y=ClaRa.Basics.Functions.Stepsmoother(func, noFunc, x)*u;

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-80,90},{-88,68},{-72,68},{-80,90}},
          lineColor={221,222,223},
          fillColor={221,222,223},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,78},{-80,-90}}, color={221,222,223}),
        Line(points={{-90,-80},{82,-80}}, color={221,222,223}),
        Polygon(
          points={{90,-80},{68,-72},{68,-88},{90,-80}},
          lineColor={221,222,223},
          fillColor={221,222,223},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{20,2},{80,-58}},
          lineColor={221,222,223},
          textString="SM"),
        Line(
          points={{-34,-82},{-34,62}},
          color={167,25,48}),
        Line(
          points={{-80,0},{-40,-40},{14,12}},
          color={221,222,223},
          smooth=Smooth.Bezier),
        Line(
          points={{-78,-80},{-34,-80},{-18,-80},{-8,-10},{14,12},{14,12},{46,34},{78,2}},
          color={27,36,42},
          smooth=Smooth.Bezier),
        Line(
          points={{14,-80},{14,64}},
          color={167,25,48})}),    Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end StepSmootherGain;
