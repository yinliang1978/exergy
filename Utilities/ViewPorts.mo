within Exergy.Utilities;
connector ViewPorts
  extends ViewPort;
  annotation (Icon(graphics={                        Ellipse(
          extent={{-72,48},{72,146}},
          lineColor={28,108,200},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),            Ellipse(
          extent={{-64,-154},{80,-56}},
          lineColor={28,108,200},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid)}), Diagram(graphics={
                                         Ellipse(
          extent={{-56,36},{50,106}},
          lineColor={28,108,200},
          fillPattern=FillPattern.Solid,
          fillColor={0,0,255}),          Ellipse(
          extent={{-50,-110},{56,-40}},
          lineColor={28,108,200},
          fillPattern=FillPattern.Solid,
          fillColor={0,0,255})}));
end ViewPorts;
