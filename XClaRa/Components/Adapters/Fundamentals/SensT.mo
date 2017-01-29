within Exergy.XClaRa.Components.Adapters.Fundamentals;
model SensT "ThermoPower's temperature sensor"

//_____________________________________________________________________________________
//  This definition was taken from Francesco Casella's library ThermoPower
//  http://sourceforge.net/projects/thermopower/
//_____________________________________________________________________________________

  replaceable package Medium = Modelica.Media.Water.WaterIF97_ph constrainedby
    Modelica.Media.Interfaces.PartialMedium "Medium model";
  Medium.BaseProperties fluid;
  FlangeA inlet(redeclare package Medium = Medium)
                                  annotation (extent=[-80, -60; -40, -20]);
  FlangeB outlet(redeclare package Medium = Medium)
                                   annotation (extent=[40, -60; 80, -20]);
  Modelica.Blocks.Interfaces.RealOutput T
    annotation (extent=[60, 40; 100, 80]);
equation
  inlet.w + outlet.w = 0 "Mass balance";
  inlet.p = outlet.p "No pressure drop";
  // Set fluid properties
  fluid.p=inlet.p;
  fluid.h = if inlet.w >= 0 then inlet.hBA else inlet.hAB;
  T = fluid.T;

  // Boundary conditions
  inlet.hAB = outlet.hAB;
  inlet.hBA = outlet.hBA;

  annotation (
    Diagram,
    Icon(Text(
        extent=[-40, 84; 38, 34],
        style(color=0, fillPattern=1),
        string="T"), graphics={
        Rectangle(
          extent={{-40,-20},{40,-60}},
          lineColor={0,0,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Ellipse(extent={{-40,100},{40,20}}, lineColor={0,0,0}),
        Line(
          points={{0,20},{0,-20}},
          color={0,0,0},
          smooth=Smooth.None),
        Line(
          points={{60,60},{40,60}},
          color={0,0,0},
          smooth=Smooth.None),
        Text(
          extent={{-40,86},{40,32}},
          lineColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="T")}),
    Documentation(info="<HTML>
<p>This component can be inserted in a hydraulic circuit to measure the temperature of the fluid flowing through it.
<p>Flow reversal is supported.
</HTML>",
      revisions="<html>
<ul>
<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Standard medium definition added.</li>
<li><i>1 Jul 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Adapted to Modelica.Media.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
end SensT;
