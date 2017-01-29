within Exergy.XClaRa.Components.Furnace.GeneralTransportPhenomena.ParticleMigration;
partial model PartialMigrationSpeed
  "Base class for the particle migration speed"

outer ClaRa.Basics.Units.VolumeFlowRate
                                    V_flow_flueGas_in;
  outer ClaRa.Basics.Units.VolumeFlowRate
                                      V_flow_flueGas_out;
  outer ClaRa.Basics.Units.Area
                            A_cross;
ClaRa.Basics.Units.Velocity
                        w;

  annotation (Icon(graphics={
        Ellipse(extent={{-100,100},{100,-100}}, lineColor={0,0,255}),
        Ellipse(
          extent={{0,30},{10,20}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-42,22},{-32,12}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-4,-12},{6,-22}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-32,-34},{-22,-44}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{26,-48},{36,-58}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{50,-16},{60,-26}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{32,4},{42,-6}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-42,58},{-32,48}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-62,-6},{-52,-16}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{-38,-64},{-28,-74}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Ellipse(
          extent={{22,76},{32,66}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,60},{-20,20},{-6,20},{-6,-60},{6,-60},{6,20},{20,20},{0,60}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-20,88},{20,58}},
          lineColor={0,0,0},
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="w")}), Documentation(info="<html>
<p><b>Model description: </b>Base class for migration speed calculation models</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));

end PartialMigrationSpeed;
