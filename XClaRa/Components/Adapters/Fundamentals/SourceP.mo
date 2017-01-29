within Exergy.XClaRa.Components.Adapters.Fundamentals;
model SourceP "Pressure source for water/steam flows"

//_____________________________________________________________________________________
//  This definition was taken from Francesco Casella's library ThermoPower
//  http://sourceforge.net/projects/thermopower/
//_____________________________________________________________________________________

  import Modelica.SIunits.*;

  type HydraulicResistance = Real (final quantity="HydraulicResistance", final unit=
             "Pa/(kg/s)");
  replaceable package Medium = Modelica.Media.Water.WaterIF97_ph constrainedby
    Modelica.Media.Interfaces.PartialMedium "Medium model";

  parameter Pressure p0=1.01325e5 "Nominal pressure";
  parameter HydraulicResistance R=0 "Hydraulic resistance";
  parameter SpecificEnthalpy h=1e5 "Nominal specific enthalpy";
  Pressure p "Actual pressure";
  FlangeB flange(redeclare package Medium=Medium)
                 annotation (extent=[80, -20; 120, 20]);
  Modelica.Blocks.Interfaces.RealInput in_p0
    annotation (extent=[-60, 72; -20, 112], rotation=-90);
  Modelica.Blocks.Interfaces.RealInput in_h
    annotation (extent=[20, 70; 60, 110], rotation=-90);
equation
  if R == 0 then
    flange.p = p;
  else
    flange.p = p + flange.w*R;
  end if;

  p = in_p0;
  if cardinality(in_p0)==0 then
    in_p0 = p0 "Pressure set by parameter";
  end if;

  flange.hBA =in_h;
  if cardinality(in_h)==0 then
    in_h = h "Enthalpy set by parameter";
  end if;
  annotation (
    Diagram,
    Icon(Text(extent=[-106, 90; -52, 50], string="p0"), Text(extent=[66, 90;
             98, 52], string="h"),
      graphics={Ellipse(
          extent={{-80,80},{80,-80}},
          lineColor={0,0,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid), Text(
          extent={{20,34},{-28,-26}},
          lineColor={255,255,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid,
          textString="P")}),
    Documentation(info="<HTML>
<p><b>Modelling options</b></p>
<p>If <tt>R</tt> is set to zero, the pressure source is ideal; otherwise, the outlet pressure decreases proportionally to the outgoing flowrate.</p>
<p>If the <tt>in_p0</tt> connector is wired, then the source pressure is given by the corresponding signal, otherwise it is fixed to <tt>p0</tt>.</p>
<p>If the <tt>in_h</tt> connector is wired, then the source pressure is given by the corresponding signal, otherwise it is fixed to <tt>h</tt>.</p>
</HTML>",
      revisions="<html>
<ul>
<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Medium model and standard medium definition added.</li>
<li><i>18 Jun 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Removed <tt>p0_fix</tt> and <tt>hfix</tt>; the connection of external signals is now detected automatically.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"));
end SourceP;
