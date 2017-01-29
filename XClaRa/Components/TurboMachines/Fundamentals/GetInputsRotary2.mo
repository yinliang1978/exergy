within Exergy.XClaRa.Components.TurboMachines.Fundamentals;
model GetInputsRotary2 "Get enabled inputs and parameters of disabled inputs"
  extends Modelica.Blocks.Interfaces.BlockIcon;

  Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft_a
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}}, rotation=
           0)));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b shaft_b
    annotation (Placement(transformation(extent={{90,-10},{110,10}},   rotation=
           0)));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}})));
end GetInputsRotary2;
