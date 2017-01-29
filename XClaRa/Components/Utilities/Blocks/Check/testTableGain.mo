within Exergy.XClaRa.Components.Utilities.Blocks.Check;
model testTableGain
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

  TableGain tableGain(table=[0,0; 1,5; 2,10])
    annotation (Placement(transformation(extent={{0,-2},{20,18}})));
  Modelica.Blocks.Sources.Ramp sine(
    startTime=2,
    height=1,
    duration=1,
    offset=1)
    annotation (Placement(transformation(extent={{-78,-2},{-58,18}})));
  Modelica.Blocks.Tables.CombiTable1D tableGain1(table=[0,0; 1,5; 2,10])
    annotation (Placement(transformation(extent={{0,-42},{20,-22}})));
  Modelica.Blocks.Math.Product product1
    annotation (Placement(transformation(extent={{30,-76},{50,-56}})));
equation
  connect(sine.y, tableGain.u[1]) annotation (Line(
      points={{-57,8},{-2,8}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(sine.y, tableGain1.u[1]) annotation (Line(
      points={{-57,8},{-30,8},{-30,-32},{-2,-32}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(tableGain1.y[1], product1.u1) annotation (Line(
      points={{21,-32},{24,-32},{24,-60},{28,-60}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(product1.u2, sine.y) annotation (Line(
      points={{28,-72},{-8,-72},{-8,-70},{-57,-70},{-57,8}},
      color={0,0,127},
      smooth=Smooth.None));
  annotation (Diagram(graphics));
end testTableGain;
