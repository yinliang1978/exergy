within Exergy.XClaRa.Components.Furnace.GeneralTransportPhenomena.BurningTime;
partial model PartialBurningTime "Base class for the burning time"

ClaRa.Basics.Units.Time
                    t;

  annotation (Icon(graphics={Ellipse(extent={{-100,100},{100,-100}}, lineColor={
              0,0,255}), Text(
          extent={{-100,10},{100,-10}},
          lineColor={0,0,255},
          textString="t_burning")}), Documentation(info="<html>
<p><b>Model description: </b>Base class for burning time models</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));
end PartialBurningTime;
