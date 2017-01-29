within Exergy.XClaRa.Components.Utilities.Blocks;
model Noise "Adds a normally distributed noise to a given mean value"
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

  parameter Real Tau_sample = 1 "Sample period" annotation(Dialog(group="Large-Scale Time Definition"));
  parameter Real startTime = 1 "start time moment" annotation(Dialog(group="Large-Scale Time Definition"));
  parameter Boolean varMeanValue=false
    "True, if mean value is time dependent input" annotation(Dialog(group="Noise Definition"));
  parameter Boolean varStandardDeviation=false
    "True, if standard deviation is time dependent input" annotation(Dialog(group="Noise Definition"));
  parameter Real mean_const= 1 "Constant mean value" annotation(Dialog(group="Noise Definition", enable=not varMeanValue));
  parameter Real stdDev_const= 1 "Constant standard deviation" annotation(Dialog(group="Noise Definition", enable=not varStandardDeviation));
  parameter Modelica.SIunits.Time Tau_smooth=Tau_sample/5
    "Time constant for smoothing of the random values" annotation(Dialog(group="Noise Definition"));

protected
  Real seed[3];
  Real y_;
  Real x;
public
  Real m_in;
  Real stdDev_in;

public
  Modelica.Blocks.Interfaces.RealOutput y annotation (Placement(transformation(extent={{100,-10},
            {120,10}}), iconTransformation(extent={{100,-10},{120,10}})));

  Modelica.Blocks.Interfaces.RealInput mean(value=m_in) if (varMeanValue)
    "Variable mass flow rate"
    annotation (Placement(transformation(extent={{-120,40},{-80,80}})));
  Modelica.Blocks.Interfaces.RealInput sigma(value=stdDev_in) if (varStandardDeviation)
    "Variable standard deviation"
    annotation (Placement(transformation(extent={{-120,-80},{-80,-40}}),
        iconTransformation(extent={{-120,-62},{-80,-22}})));
algorithm

equation
  if (not varMeanValue) then
    m_in=mean_const;
  end if;
  if (not varStandardDeviation) then
    stdDev_in=stdDev_const;
  end if;

   x=smooth(1,div(time,Tau_sample));

  seed =  {noEvent(smooth(1,x)), noEvent(smooth(1,x)^3), noEvent(smooth(1,x^2))};
  if noEvent(time > startTime) then
    der(y_) = (
      Exergy.XClaRa.Components.Utilities.Blocks.Fundamentals.normalvariate(
                m_in,
                stdDev_in,
                seed) - y_)/Tau_smooth;
  else
    der(y_)=(m_in-y_)/Tau_smooth;
  end if;
  der(y)=(y_-y)/Tau_smooth;
initial equation
  y=m_in;
  y_=m_in;
  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(points={{-78,82},{-78,-76},{78,-76}}, color={221,222,223}),
        Line(
          points={{-78,0},{-60,-20},{-42,30},{-32,-14},{-20,-14},{-4,-64},{4,24},{8,-28},{24,76},{36,-16},{38,-54},{69.9688,33.9141},{72,-12}},
          color={27,36,42},
          smooth=Smooth.Bezier)}),
    Documentation(info="<html>
This random generator combines ideas from:
<p>
<b>Peter Fritzson</b>: \"Principles of Object-Oriented Modeling and Simulation with Modelica 2.1\" published by Wiley Interscience, 2004. <br>
<b>Dorel Aiordachioaie et al.</b>: \"On Noise Modelling and Simulation\" In Proceedings of the Modelica Conference 2006.
</p>
</html>", revisions="<html>

<ul>
<li> version 0.1: <b>Friedrich Gottelt, XRG Simulation GmbH</b> Initial implimentation. 
</ul>
</html>"));
end Noise;
