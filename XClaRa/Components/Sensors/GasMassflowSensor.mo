within Exergy.XClaRa.Components.Sensors;
model GasMassflowSensor "Ideal one port mass flow sensor"
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

extends Exergy.XClaRa.Components.Sensors.GasSensorBase;
  outer ClaRa.SimCenter simCenter;

  parameter Integer unitOption = 1 "Unit of output" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"), choices(choice=1 "kg/s", choice=2 "t/h", choice=3 "kg/h", choice=4 "t/s"));
  Modelica.Blocks.Interfaces.RealOutput m_flow(final quantity="mass flow",
                                                            displayUnit = "kg/s",
    final unit="kg/s") "mass flow in port"
    annotation (Placement(transformation(extent={{100,-10},{120,10}},
                                                                    rotation=
            0), iconTransformation(extent={{100,-10},{120,10}})));

equation
  if unitOption==1 then // kg/s
    m_flow = inlet.m_flow;
  elseif unitOption==2 then // tons per hour
    m_flow = inlet.m_flow*3.6;
  elseif unitOption==3 then // kg/h
    m_flow = inlet.m_flow*3600;
  elseif unitOption==4 then // tons per s
    m_flow = inlet.m_flow/1000;
  else
    m_flow=-1; //dummy
    assert(false, "Unknown unit option in " + getInstanceName());
  end if;

  annotation (Diagram(graphics), Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                                      graphics={
        Text(
          extent={{-100,40},{100,0}},
          lineColor={27,36,42},
          fillColor={0,255,0},
          fillPattern=FillPattern.Solid,
          textString="FIT"),
        Text(
          extent={{-100,0},{100,-40}},
          lineColor={27,36,42},
          fillColor={0,255,0},
          fillPattern=FillPattern.Solid,
          textString="%name"),
        Text(
          extent={{-100,60},{100,90}},
          fillColor={215,215,215},
          fillPattern=FillPattern.Solid,
          lineColor=DynamicSelect({230, 230, 230},  if p > 0 then {0,131,169} else {167,25,48}),
          textString=DynamicSelect(" m_flow ", String(m_flow,significantDigits=integer(1))+" kg/s")),
        Line(
          points={{-98,-100},{96,-100}},
          color={118,106,98},
          smooth=Smooth.None,
          thickness=0.5)}));
end GasMassflowSensor;
