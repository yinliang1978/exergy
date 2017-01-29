within Exergy.XClaRa.Components.Utilities.Blocks;
model Convert2StateDiscrete
  "Converts a flaoting value to a discrete one. Value is changed when a certain threshold is violated"
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
  parameter Real threshold = 0.1 "Threshold of identification of changes";
  parameter Boolean smoothValue = true
    "True, if output value should be smoothed using a firstOrder";
  parameter ClaRa.Basics.Units.Time Tau=0.1 "Time constant for smoothing" annotation(Dialog(enable=smoothValue));
  discrete Real y_;
  discrete Real py_;
equation

  when abs(pre(y_)-u)>threshold then
    y_=u;
    py_=pre(y_);
  end when;

  if smoothValue then
    der(y) = (y_-y)/Tau;
  else
    y=y_;
  end if;
initial equation
  y_=u;
  if smoothValue then
    y=u;
  end if;

equation
  connect(y, y) annotation (Line(
      points={{110,0},{106,0},{106,0},{110,0}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
                             Line(
          points={{-100,0},{-50,60},{30,-60},{100,20}},
          color={221,222,223},
          smooth=Smooth.Bezier), Line(
          points={{-100,0},{-80,0},{-80,22},{-12,22},{-12,0},{10,0},{10,-26},{
              82,-26},{82,0},{100,0}},
          color={27,36,42},
          smooth=Smooth.None)}));
end Convert2StateDiscrete;
