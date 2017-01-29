within Exergy.XThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor;
model Mechanical_losses
  "Computes the Mechanical losses as a sum of different contributions (constant losses, proportional losses, friction torque)"
parameter Modelica.SIunits.Power W_dot_loss_0 = 0 "Constant losses";
parameter Real alpha = 0 "Proportionality factor (Wdot_loss = alpha*Wdot_tot)";
parameter Modelica.SIunits.Torque T_friction= 0
    "Friction torque (Wdot_loss = 2*pi*N_rot*T_friction)";
Modelica.SIunits.Power W_dot_loss;
Modelica.SIunits.Power W_dot_propor;
Modelica.SIunits.Power W_dot_friction;
Modelica.SIunits.Power W_dot_A;
Modelica.SIunits.Power W_dot_B;
Modelica.SIunits.AngularVelocity omega "Angular velocity of the shaft [rad/s] ";

  Modelica.Mechanics.Rotational.Interfaces.Flange_a flange_A      annotation (
      Placement(transformation(extent={{-116,-16},{-84,16}}),
                                                          iconTransformation(
          extent={{-112,-10},{-92,10}})));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_B      annotation (
      Placement(transformation(extent={{84,-16},{116,16}}),
                                                          iconTransformation(
          extent={{90,-10},{110,10}})));
  Modelica.Blocks.Interfaces.RealOutput Wdot_loss  annotation (Placement(transformation(
        extent={{-16,-16},{16,16}},
        rotation=0,
        origin={98,-68}),iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={82,-50})));
equation
omega = der(flange_A.phi);

W_dot_A = omega  * flange_A.tau;
W_dot_B = -omega  * flange_B.tau;
flange_A.phi = flange_B.phi;

W_dot_loss = W_dot_loss_0 + W_dot_propor + W_dot_friction;
Wdot_loss = W_dot_loss;

W_dot_propor = alpha* abs(W_dot_A);
W_dot_friction = abs(omega*T_friction);

W_dot_B = W_dot_A - W_dot_loss;

  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics={
        Rectangle(
          extent={{-102,10},{98,-10}},
          lineColor={0,0,0},
          fillPattern=FillPattern.HorizontalCylinder,
          fillColor={192,192,192}),
        Rectangle(extent={{-62,-10},{58,-60}}, lineColor={0,0,0}),
        Rectangle(
          extent={{-62,-10},{58,-25}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-62,-45},{58,-61}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-52,-18},{48,-50}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{58,-60},{58,-70},{73,-70},{73,-80},{-77,-80},{-77,-70},{-62,-70},
              {-62,-60},{58,-60}},
          lineColor={0,0,0},
          fillColor={160,160,164},
          fillPattern=FillPattern.Solid),
        Line(points={{-77,-10},{-77,-70}}, color={0,0,0}),
        Line(points={{73,-10},{73,-70}}, color={0,0,0}),
        Rectangle(extent={{-62,60},{58,10}}, lineColor={0,0,0}),
        Rectangle(
          extent={{-62,60},{58,45}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-62,25},{58,10}},
          lineColor={0,0,0},
          fillColor={192,192,192},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-52,51},{48,19}},
          lineColor={0,0,0},
          fillColor={255,255,255},
          fillPattern=FillPattern.Solid),
        Line(points={{-77,70},{-77,10}}, color={0,0,0}),
        Polygon(
          points={{58,60},{58,70},{73,70},{73,80},{-77,80},{-77,70},{-62,70},{-62,
              60},{58,60}},
          lineColor={0,0,0},
          fillColor={160,160,164},
          fillPattern=FillPattern.Solid),
        Line(points={{73,70},{73,10}}, color={0,0,0}),
        Text(
          extent={{-152,130},{148,90}},
          textString="%name",
          lineColor={0,0,255})}));
end Mechanical_losses;
