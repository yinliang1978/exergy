within Exergy.Utilities;
model ViewRoute
  parameter Integer nSubSystem=0
   annotation(Evaluate=true, Dialog(connectorSizing=true, tab="General",group="Ports"));

  Utilities.ViewPorts viewOutput[nSubSystem]
    annotation (Placement(transformation(extent={{-128,-32},{-74,28}})));
  Utilities.ViewPort viewTotal
    annotation (Placement(transformation(extent={{76,-28},{120,26}})));

equation
  viewTotal.E.E = sum(viewOutput.E.E);
  viewTotal.E.Ex = sum(viewOutput.E.Ex);
  //  viewTotal.internalEnergy.E_int = sum(viewOutput.internalEnergy.E_int);
  //viewTotal.internalEnergy.Ex_int = sum(viewOutput.internalEnergy.Ex_int);

 //   viewTotal.internalEnergy.E = der(viewTotal.internalEnergy.E_int);
 // viewTotal.internalEnergy.Ex = der(viewTotal.internalEnergy.Ex_int);
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{
            -100,-100},{100,100}}), graphics={Rectangle(
          extent={{-70,50},{40,-58}},
          lineColor={28,108,200},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid), Polygon(
          points={{40,32},{78,0},{40,-40},{40,32}},
          lineColor={28,108,200},
          fillColor={170,255,255},
          fillPattern=FillPattern.Solid)}), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
        graphics={Rectangle(
          extent={{-68,50},{42,-58}},
          lineColor={28,108,200},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid), Polygon(
          points={{42,32},{80,0},{42,-40},{42,32}},
          lineColor={28,108,200},
          fillColor={170,255,255},
          fillPattern=FillPattern.Solid)}));
end ViewRoute;
