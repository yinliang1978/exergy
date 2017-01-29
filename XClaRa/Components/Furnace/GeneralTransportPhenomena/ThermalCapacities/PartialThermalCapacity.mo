within Exergy.XClaRa.Components.Furnace.GeneralTransportPhenomena.ThermalCapacities;
partial model PartialThermalCapacity
  import ClaRa;
  ClaRa.Basics.Interfaces.HeatPort_a   heat_in
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat_out
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

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
<p><b>Model description: </b>Base class for thermal capacities</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));

end PartialThermalCapacity;
