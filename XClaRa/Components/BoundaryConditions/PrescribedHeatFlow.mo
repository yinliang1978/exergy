within Exergy.XClaRa.Components.BoundaryConditions;
model PrescribedHeatFlow "Prescribed heat flow boundary condition 1D"
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                        //
//                                                                           //
// Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
// Copyright © 2013-2016, DYNCAP/DYNSTART research team.                     //
//___________________________________________________________________________//
// DYNCAP and DYNSTART are research projects supported by the German Federal //
// Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
// The research team consists of the following project partners:             //
// Institute of Energy Systems (Hamburg University of Technology),           //
// Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
// TLK-Thermo GmbH (Braunschweig, Germany),                                  //
// XRG Simulation GmbH (Hamburg, Germany).                                   //
//___________________________________________________________________________//

  parameter Modelica.SIunits.Length length "Length of cylinder" annotation(Dialog(group="Geometry"));
  parameter Integer N_axial = 3 "Number of axial elements" annotation(Dialog(group="Discretisation"));
  parameter Modelica.SIunits.Length Delta_x[N_axial]=ClaRa.Basics.Functions.GenerateGrid(        {1,-1}, length, N_axial)
    "Discretisation scheme"  annotation(Dialog(group="Discretisation"));

  Modelica.Blocks.Interfaces.RealInput Q_flow
        annotation (Placement(transformation(
        origin={-100,0},
        extent={{20,-20},{-20,20}},
        rotation=180)));
  ClaRa.Basics.Interfaces.HeatPort_b port[N_axial]
                             annotation (Placement(transformation(extent={{90,
            -10},{110,10}}, rotation=0)));
equation
  port.Q_flow = -Q_flow.*(Delta_x/sum(Delta_x));
  annotation (
    Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}}), graphics={
        Line(
          points={{-60,-20},{40,-20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-60,20},{40,20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-80,0},{-60,-20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-80,0},{-60,20}},
          color={191,0,0},
          thickness=0.5),
        Polygon(
          points={{40,0},{40,40},{70,20},{40,0}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{40,-40},{40,0},{70,-20},{40,-40}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{70,40},{90,-40}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid)}),
    Documentation(info="<HTML>
<p>
This model allows a specified amount of heat flow rate to be \"injected\"
into a thermal system at a given port.  The amount of heat
is given by the input signal Q_flow into the model. The heat flows into the
component to which the component PrescribedHeatFlow is connected,
if the input signal is positive.
</p>
<p>
If parameter alpha is > 0, the heat flow is mulitplied by (1 + alpha*(port.T - T_ref))
in order to simulate temperature dependent losses (which are given an reference temperature T_ref).
</p>
</HTML>
"), Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,
            100}}),      graphics={
        Line(
          points={{-60,-20},{68,-20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-60,20},{68,20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-80,0},{-60,-20}},
          color={191,0,0},
          thickness=0.5),
        Line(
          points={{-80,0},{-60,20}},
          color={191,0,0},
          thickness=0.5),
        Polygon(
          points={{60,0},{60,40},{90,20},{60,0}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{60,-40},{60,0},{90,-20},{60,-40}},
          lineColor={191,0,0},
          fillColor={191,0,0},
          fillPattern=FillPattern.Solid)}));
end PrescribedHeatFlow;
