within Exergy.XThermoCycle.Components.Units.Tanks;
package HeatStorageWaterHeater "Stratified tank with an internal heat exchanger and ambient heat losses (nodal model)"





  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics={
        Ellipse(
          extent={{-42,52},{38,80}},
          lineColor={0,0,0},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-42,64},{38,-92}},
          lineColor={215,215,215},
          fillPattern=FillPattern.Solid,
          fillColor={215,215,215}),
        Line(
          points={{38,26},{-12,26},{18,14},{-12,4},{18,-6},{-10,-18},{20,-28},{
              -10,-40},{36,-40}},
          color={0,0,0},
          smooth=Smooth.None,
          thickness=0.5),
        Line(
          points={{26,-44},{38,-40},{26,-36}},
          color={0,0,0},
          thickness=0.5,
          smooth=Smooth.None),
        Line(
          points={{-42,-92},{38,-92},{38,66}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-42,66},{-42,-92}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{-42,58},{38,58}},
          color={0,0,0},
          smooth=Smooth.None)}), Diagram(coordinateSystem(preserveAspectRatio=
            false, extent={{-100,-100},{100,100}}), graphics),
    Documentation(info="<html>
<p>Nodal model of a stratified tank, with the following hypotheses:</p>
<p><ul>
<li>No heat transfer between the different nodes</li>
<li>The internal heat exchanger is discretized in the same way as the tank: each cell of the heat exchanger corresponds to one cell of the tank and exchanges heat with that cell only.</li>
<li>Incompressible fluid in both the tank and the heat exchanger</li>
</ul></p>
<p><br/>The tank is discretized using a modified version of the incompressible Cell1Dim model adding an additional heat port. The heat exchanger is modeled using the Flow1Dim component and a wall component.</p>
</html>"));
end HeatStorageWaterHeater;
