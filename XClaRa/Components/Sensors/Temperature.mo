within Exergy.XClaRa.Components.Sensors;
model Temperature "Ideal one port temperature sensor"
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

  extends ClaRa.Basics.Icons.Sensor1;
  outer ClaRa.SimCenter simCenter;
  parameter TILMedia.VLEFluidTypes.BaseVLEFluid   medium= simCenter.fluid1
    "Medium to be used"                                                                        annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));
  parameter Integer unitOption = 1 "Unit of output" annotation(choicesAllMatching, Dialog( group="Fundamental Definitions"), choices(choice=1 "Kelvin", choice=2
        "Degree Celsius",                                                                     choice=3
        "Degree Fahrenheit"));
  ClaRa.Basics.Units.Temperature_DegC T_celsius "Temperatur in Degree Celsius";
  Modelica.Blocks.Interfaces.RealOutput T "Temperature in port medium"
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
                                                                    rotation=
            0), iconTransformation(extent={{100,-10},{120,10}})));

public
  ClaRa.Basics.Interfaces.FluidPortIn port(Medium=medium)
    annotation (Placement(transformation(extent={{-10,-110},{10,-90}}),
        iconTransformation(extent={{-10,-110},{10,-90}})));
protected
  TILMedia.VLEFluid_ph fluid(
    h=inStream(port.h_outflow),
    p=port.p,
    xi=inStream(port.xi_outflow),
    vleFluidType=medium)
    annotation (Placement(transformation(extent={{-48,28},{-28,48}})));
equation
  if unitOption == 1 then //Kelvin
    T = fluid.T;
  elseif unitOption == 2 then // Degree Celsius
    T = Modelica.SIunits.Conversions.to_degC(fluid.T);
  elseif unitOption == 3 then // Degree Fahrenheit
    T = Modelica.SIunits.Conversions.to_degF(fluid.T);
  else
    T=-1;  //dummy
    assert(false, "Unknown unit option in " + getInstanceName());
  end if;

  T_celsius = fluid.T - 273.15;

  port.m_flow = 0;
  port.h_outflow = 0;
  port.xi_outflow = zeros(medium.nc-1);
 // port.C_outflow = zeros(Medium.nC);

  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                                      graphics={
        Text(
          extent={{-100,0},{100,-40}},
          lineColor={27,36,42},
          fillColor={0,255,0},
          fillPattern=FillPattern.Solid,
          textString="%name"),
        Text(
          extent={{-100,40},{100,0}},
          lineColor={27,36,42},
          fillColor={0,255,0},
          fillPattern=FillPattern.Solid,
          textString="TIT"),
        Text(
          extent={{-100,60},{60,90}},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          lineColor=DynamicSelect({230, 230, 230},  if T_celsius > 0 then {0,131,169} else {167,25,48}),
          textString=DynamicSelect(" T ", String(T, format="1.1f"))),
        Text(
          extent={{50,90},{90,60}},
          lineColor=DynamicSelect({230, 230, 230},  if T_celsius>0 then {0,131,169} else {167,25,48}),
          textString=DynamicSelect("", if unitOption==1 then "K" elseif unitOption==2 then "°C" else "°F"),
          horizontalAlignment=TextAlignment.Left)}));
end Temperature;
