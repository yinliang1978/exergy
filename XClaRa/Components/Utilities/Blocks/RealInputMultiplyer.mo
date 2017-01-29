within Exergy.XClaRa.Components.Utilities.Blocks;
block RealInputMultiplyer
  "Distributes a single real input into N real outputs having the same value"
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

  parameter Integer N=1 "Number outputs";

  Modelica.Blocks.Interfaces.RealInput Signal "Input signal to be duplicated"
    annotation (Placement(transformation(extent={{-136,-18},{-100,18}}),
        iconTransformation(extent={{-136,-18},{-100,18}})));

  Modelica.Blocks.Interfaces.RealOutput y[N]
    annotation (Placement(transformation(extent={{100,-11},{120,10}}), iconTransformation(extent={{100,-11},{120,10}})));

equation
  y = {Signal for i in 1:N};

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                   graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-98,-40},{-18,-84}},
          lineColor={221,222,223},
          textString="N = %N"),
        Ellipse(
          extent={{-14,16},{16,-14}},
          lineColor={221,222,223},
          fillColor={221,222,223},
          fillPattern=FillPattern.Solid),
        Line(points={{-100,0},{-6,0}}, color={221,222,223}),
        Line(points={{0,0},{100,10}}, color={221,222,223}),
        Line(points={{100,0},{10,0}}, color={221,222,223}),
        Line(points={{0,0},{100,-10}}, color={221,222,223})}),
                                 Diagram(graphics));
end RealInputMultiplyer;
