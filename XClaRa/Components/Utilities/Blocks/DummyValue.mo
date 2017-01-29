within Exergy.XClaRa.Components.Utilities.Blocks;
model DummyValue "A dummy connection for Real inputs"
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

  Modelica.Blocks.Interfaces.RealOutput y
    annotation (Placement(transformation(extent={{100,-10},{120,10}}), iconTransformation(extent={{100,-10},{120,10}})));

equation
  y=0;
  annotation (Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        initialScale=0.01), graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          pattern=LinePattern.Solid,
          fillColor={27,36,42},
          fillPattern=FillPattern.Solid,
          lineColor={27,36,42},
          radius=20),
        Rectangle(
          extent={{-100,100},{0,0}},
          pattern=LinePattern.Solid,
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid,
          lineColor={27,36,42}),
        Rectangle(
          extent={{0,0},{100,-100}},
          pattern=LinePattern.Solid,
          fillColor={235,183,0},
          fillPattern=FillPattern.Solid,
          lineColor={27,36,42})}),
                                Diagram(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100,-100},{100,100}},
        initialScale=0.01)));
end DummyValue;
