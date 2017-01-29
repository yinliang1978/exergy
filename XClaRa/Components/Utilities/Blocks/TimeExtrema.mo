within Exergy.XClaRa.Components.Utilities.Blocks;
model TimeExtrema
  "Calculates the minimum and maximum value in a given period of time"
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
  parameter SI.Time startTime= 0 "Start time for min/max evaluation";
  parameter Integer initOption= 0 "Init option |initial u| y_min/max_start" annotation(Dialog(group="Initialisation"), choices(choice=0
        "initial u",                                                                                                    choice=1
        "initial y_min/y_max"));
  parameter Real y_start[2]= {10,-10} "Y_min_start | y_max_start"  annotation(Dialog(group="Initialisation", enable=initOption==1));

protected
  Real ymax;
  Real ymin;

public
  Modelica.Blocks.Interfaces.RealOutput y_max
    annotation (Placement(transformation(extent={{100,50},{120,70}})));
  Modelica.Blocks.Interfaces.RealOutput y_min
    annotation (Placement(transformation(extent={{100,-50},{120,-30}})));
  Modelica.Blocks.Interfaces.RealInput u
    annotation (Placement(transformation(extent={{-140,-20},{-100,20}})));
equation
  der(ymax)= if  time >= startTime then noEvent(if ymax<u then (u-ymax)/1 else 0) else 0;
  der(ymin)= if  time >= startTime then noEvent(if ymin>u then (u-ymin)/1 else 0) else 0;
  y_max=ymax;
  y_min=ymin;
initial equation
  if initOption==0 then
    ymax=u;
    ymin=u;
  elseif initOption ==1 then
    ymin=y_start[1];
    ymax=y_start[2];
  else
    assert(false, "Unknown initi option in component " + getInstanceName());
  end if;
  annotation (Icon(graphics={
        Rectangle(
          extent={{-100,102},{100,-98}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-80,80},{-80,-78},{76,-78}},
          color={221,222,223},
          smooth=Smooth.None),
        Line(
          points={{-80,-16},{-56,34},{-30,64},{-26,-36},{-6,-38},{-2,-28},{16,20},{43.9805,-7.98047},{76,2}},
          color={221,222,223},
          smooth=Smooth.Bezier),
        Line(
          points={{-80,-18},{-58,-18},{-28,-18},{-28,-18},{-22,-40},{4,-40},{4,-40},{76,-40}},
          color={27,36,42},
          smooth=Smooth.Bezier),
        Line(
          points={{-80,-14},{-56,36},{-38,58},{-4,58},{76,58}},
          color={167,25,48},
          smooth=Smooth.Bezier)}),
                                 Diagram(graphics));
end TimeExtrema;
