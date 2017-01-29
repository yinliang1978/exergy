within Exergy.XClaRa.Components.Utilities.Blocks;
block VarGain
  "Output the product of a variable gain value with the input signal"
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
  input Real k(start=1, unit="1")
    "Variable gain value multiplied with input signal"                               annotation(Dialog);
public
  Modelica.Blocks.Interfaces.RealInput u "Input signal connector"
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}},
      rotation=0)));
  Modelica.Blocks.Interfaces.RealOutput y "Output signal connector"
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
      rotation=0)));

equation
  y = k*u;
  annotation (
    Documentation(info="
<HTML>
<p>
This block computes output <i>y</i> as
<i>product</i> of gain <i>k</i> with the
input <i>u</i>:
</p>
<pre>
    y = k * u;
</pre>

</HTML>
"), Icon(coordinateSystem(
    preserveAspectRatio=true,
    extent={{-100,-100},{100,100}},
    grid={2,2}), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
    Polygon(
          points={{-100,-100},{-100,100},{100,0},{-100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
    Text(
      extent={{-150,-140},{150,-100}},
      lineColor={0,0,0},
      textString="k=%k"),
    Text(
      extent={{-150,140},{150,100}},
      textString="%name",
      lineColor={0,48,111}),
        Line(
          points={{-100,-60},{60,20}},
          color={221,222,223})}),
    Diagram(coordinateSystem(
    preserveAspectRatio=true,
    extent={{-100,-100},{100,100}},
    grid={2,2}), graphics={Polygon(
      points={{-100,-100},{-100,100},{100,0},{-100,-100}},
      lineColor={0,0,127},
      fillColor={255,255,255},
      fillPattern=FillPattern.Solid), Text(
      extent={{-76,38},{0,-34}},
      textString="k",
      lineColor={0,0,255})}));
end VarGain;
