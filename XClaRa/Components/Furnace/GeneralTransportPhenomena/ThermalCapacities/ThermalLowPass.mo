within Exergy.XClaRa.Components.Furnace.GeneralTransportPhenomena.ThermalCapacities;
model ThermalLowPass "A simple thermal low pass with time constant tau"
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

  extends
    Exergy.XClaRa.Components.Furnace.GeneralTransportPhenomena.ThermalCapacities.PartialThermalCapacity;

public
  parameter ClaRa.Basics.Units.Time Tau=1 "Time constant for thermal low pass";
  parameter ClaRa.Basics.Units.Temperature T_out_initial=273.15 + 250
    "Initial value for Temperature at heat_out";

initial equation
  heat_out.T  = T_out_initial;

equation
  der(heat_out.T) = 1/Tau*(heat_in.T - heat_out.T);
  heat_out.Q_flow  + heat_in.Q_flow = 0;
  annotation (Diagram(graphics), Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-80,92},{-80,-72}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-82,-70},{82,-70}},
          color={0,0,0},
          smooth=Smooth.None),
        Polygon(
          points={{77,-70},{77,-72},{77,-68},{77,-72},{83,-70},{77,-68},{77,-70}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{-3,0},{-3,-2},{-3,2},{-3,-2},{3,0},{-3,2},{-3,0}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid,
          origin={-80,90},
          rotation=90),
        Line(
          points={{-80,60},{-22,60},{78,-44}},
          color={0,0,0},
          smooth=Smooth.None)}),
    Documentation(info="<html>
<p><b>Model description: </b>Model for a thermal low pass</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));
end ThermalLowPass;
