within Exergy.XClaRa.Components.Adapters.Fundamentals;
connector FlangeA "A-type flange connector for water/steam flows"

//_____________________________________________________________________________________
//  This connector definition was taken from Francesco Casella's library ThermoPower
//  http://sourceforge.net/projects/thermopower/
//_____________________________________________________________________________________

  replaceable package Medium = Modelica.Media.Water.WaterIF97_ph constrainedby
    Modelica.Media.Interfaces.PartialMedium "Medium model";
  Medium.AbsolutePressure p "Pressure";
  flow Medium.MassFlowRate w "Mass flowrate";
  output Medium.SpecificEnthalpy hAB "Specific enthalpy of fluid going out";
  input Medium.SpecificEnthalpy hBA "Specific enthalpy of entering fluid";
  annotation (Documentation(info="<HTML>
<p> Must always be connected to a single type-B connector <tt>FlangeB</tt>.
</HTML>",
      revisions="<html>
<ul>
<li><i>16 Dec 2004</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       Medium model added.</li>
<li><i>1 Oct 2003</i>
    by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
       First release.</li>
</ul>
</html>"), Icon(graphics={Ellipse(
          extent={{-100,100},{100,-100}},
          lineColor={0,0,255},
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid)}));
end FlangeA;
