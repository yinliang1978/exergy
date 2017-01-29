within Exergy.XClaRa.Components.Utilities.Blocks;
model TrackThreshold
  "A state based threshold observation with one possible switch"
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
  //extends Modelica.Blocks.Interfaces.partialBooleanBlockIcon;
  parameter Real threshold= 1 "Threshold to be observed";
  parameter Boolean switchWhenGreater=true
    "True if output switschen when u>thershold";
  parameter Boolean y_start= false "Initial output value";

  Boolean hasSwitched(start=false);
  Modelica.Blocks.Interfaces.RealInput u annotation (Placement(transformation(
          extent={{-140,-20},{-100,20}}),iconTransformation(extent={{-140,-20},{
            -100,20}})));
  Modelica.Blocks.Interfaces.BooleanOutput y(start=y_start)
    annotation (Placement(transformation(extent={{100,-10},{120,10}})));

algorithm
  if if switchWhenGreater then u>threshold else u<threshold then
    if not hasSwitched then
      hasSwitched :=true;
      y:=not pre(y);
    else
      y:=pre(y);
    end if;

  end if;

  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-78,-50},{-48,-38},{14,34},{62,-40},{80,6}},
          color={221,222,223},
          smooth=Smooth.Bezier),
        Line(
          points={{-82,-20},{80,-20}},
          color={27,36,42},
          smooth=Smooth.Bezier),
        Line(
          points={{-82,0},{-34,0},{-34,62},{80,62}},
          color={255,0,255},
          smooth=Smooth.None)}));
end TrackThreshold;
