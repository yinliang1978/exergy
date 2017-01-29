within Exergy.XClaRa.Components.Adapters;
model Fluid2ClaRa "A Modelica.Fluid to ClaRa connector"
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

  replaceable package ModelicaMedium =
      Modelica.Media.Water.WaterIF97_ph "Medium model of the MSL part"
                                   annotation (choicesAllMatching=true, Dialog(group="Fundamental Definitions"));
  outer ClaRa.SimCenter simCenter;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid ClaRaMedium = simCenter.fluid1
    "Medium in the component"
    annotation(choices(choice=simCenter.fluid1
        "First fluid defined in global simCenter",
                       choice=simCenter.fluid2
        "Second fluid defined in global simCenter",
                       choice=simCenter.fluid3
        "Third fluid defined in global simCenter"),       Dialog(group="Fundamental Definitions"));
  Modelica.Fluid.Interfaces.FluidPort_a port_a(redeclare package Medium =
        ModelicaMedium)
    annotation (Placement(transformation(extent={{-108,-12},{-88,8}}),
        iconTransformation(extent={{-108,-10},{-88,10}})));
  ClaRa.Basics.Interfaces.FluidPortIn steam_a(Medium=ClaRaMedium)                 annotation (Placement(
        transformation(extent={{90,-9.9},{109,9.9}}), iconTransformation(extent=
           {{90,-10},{109,9.95}})));

equation
  port_a.m_flow+steam_a.m_flow = 0;
  port_a.h_outflow = inStream(steam_a.h_outflow);
  inStream(port_a.h_outflow) = steam_a.h_outflow;

  port_a.p = steam_a.p;
  inStream(port_a.Xi_outflow) =steam_a.xi_outflow;
  port_a.Xi_outflow = inStream(steam_a.xi_outflow);

  port_a.C_outflow=zeros(ModelicaMedium.nC);
  annotation (Icon(graphics={Polygon(
          points={{-96,10},{2,10},{36,-10},{-98,-10},{-96,10}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={0,128,255},
          fillPattern=FillPattern.Solid), Polygon(
          points={{100,10},{100,10},{-20,10},{-20,10},{-12,0},{0,0},{12,0},{20,
              -10},{20,-10},{100,-10},{100,-10},{100,10}},
          smooth=Smooth.Bezier,
          fillColor={0,131,169},
          fillPattern=FillPattern.Solid,
          lineColor={0,131,169})}),         Diagram(graphics));
end Fluid2ClaRa;
