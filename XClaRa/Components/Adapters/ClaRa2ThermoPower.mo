within Exergy.XClaRa.Components.Adapters;
model ClaRa2ThermoPower
  "Adapter for ClaRa to ThermoPower (v2.2) fluid connector"
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                            //
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
  replaceable package ThermoPowerMedium =
      Modelica.Media.Water.WaterIF97_ph "Medium model of the ThermoPower part"
                                           annotation (choicesAllMatching=true, Dialog(group="Fundamental Definitions"));
  outer ClaRa.SimCenter simCenter;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid ClaRaMedium = simCenter.fluid1
    "Medium model of the ClaRa part"
    annotation(choices(choice=simCenter.fluid1
        "First fluid defined in global simCenter",
                       choice=simCenter.fluid2
        "Second fluid defined in global simCenter",
                       choice=simCenter.fluid3
        "Third fluid defined in global simCenter"),       Dialog(group="Fundamental Definitions"));

  Fundamentals.FlangeB flangeB(redeclare package Medium =
        ThermoPowerMedium)                                                   annotation (Placement(transformation(extent={{-108,-12},{-88,8}}), iconTransformation(extent={{90,-10},{110,10}})));
  ClaRa.Basics.Interfaces.FluidPortIn inlet(Medium=ClaRaMedium)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));

equation
  flangeB.w+inlet.m_flow = 0;
  flangeB.hAB = inlet.h_outflow;
  flangeB.hBA = inStream(inlet.h_outflow);

  flangeB.p = inlet.p;
  inlet.xi_outflow =ones(ClaRaMedium.nc-1);

  annotation (Icon(graphics={Polygon(
          points={{98,10},{0,10},{-34,-10},{100,-10},{98,10}},
          lineColor={0,0,255},
          smooth=Smooth.None,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid), Polygon(
          points={{-100,10},{-100,10},{20,10},{20,10},{12,0},{0,0},{-12,0},{-20,
              -10},{-20,-10},{-100,-10},{-100,-10},{-100,10}},
          lineColor={0,0,0},
          smooth=Smooth.Bezier,
          fillColor={0,131,169},
          fillPattern=FillPattern.Solid)}), Diagram(graphics));
end ClaRa2ThermoPower;
