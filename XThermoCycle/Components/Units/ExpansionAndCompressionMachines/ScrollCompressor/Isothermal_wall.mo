within Exergy.XThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor;
model Isothermal_wall
Modelica.SIunits.TemperatureSlope der_T(start=0);
parameter Modelica.SIunits.HeatCapacity C = 400;
parameter Modelica.SIunits.ThermalConductance AU_amb = 10;
Modelica.SIunits.Power Q_dot_amb;
Modelica.SIunits.Temperature T_wall(start = 290);
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a Compresseur
    annotation (Placement(transformation(extent={{-16,58},{16,90}}),
        iconTransformation(extent={{-8,16},{6,30}})));
  //Pertes venant du compresseur
  Modelica.Blocks.Interfaces.RealInput W_dot_loss annotation (Placement(transformation(
        extent={{20,20},{-20,-20}},
        rotation=90,
        origin={-74,32}), iconTransformation(
        extent={{6,6},{-6,-6}},
        rotation=90,
        origin={-60,18})));
  //Température ambiante en entrée
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a Ambient
    annotation (Placement(transformation(extent={{-16,-98},{16,-66}}),
        iconTransformation(extent={{-8,-30},{6,-16}})));
equation
Q_dot_amb = AU_amb*(Ambient.T - T_wall);
Ambient.Q_flow = Q_dot_amb;
C*der_T = Compresseur.Q_flow + Q_dot_amb + W_dot_loss;
Compresseur.T = T_wall;
der_T = der(T_wall);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={      Rectangle(
          extent={{-100,16},{100,-16}},
          lineColor={0,0,0},
          fillColor={128,128,128},
          fillPattern=FillPattern.Solid), Text(
          extent={{-92,8},{98,-10}},
          lineColor={0,0,0},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          textString="Fictitious isothermal wall")}), Diagram(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end Isothermal_wall;
