within Exergy.XClaRa.Components.TurboMachines.Pumps;
partial model Pump_Base "Base class for pumps"
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
  import SI = ClaRa.Basics.Units;
  extends ClaRa.Basics.Icons.Pump;
  //extends Modelica.Icons.UnderConstruction;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid   medium=simCenter.fluid1
    "Medium in the component"            annotation(choicesAllMatching=true, Dialog(group="Fundamental Definitions"));
  parameter Boolean showExpertSummary = simCenter.showExpertSummary
    "True, if expert summary should be applied"                                             annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                            annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));
  SI.Pressure Delta_p
    "Pressure difference between pressure side and suction side";
  SI.VolumeFlowRate V_flow "Volume flow rate";
  SI.Power P_hyd;

  outer ClaRa.SimCenter simCenter;

  ClaRa.Basics.Interfaces.FluidPortIn inlet(
                                           Medium=medium)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.FluidPortOut
                                     outlet(Medium=medium)
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
protected
  TILMedia.VLEFluid_ph  fluidIn(vleFluidType =    medium, p=inlet.p,
    h=homotopy(actualStream(inlet.h_outflow), inStream(inlet.h_outflow)))
    annotation (Placement(transformation(extent={{-88,-12},{-68,8}})));
  TILMedia.VLEFluid_ph  fluidOut(
                                vleFluidType =    medium,
    p=outlet.p,
    h=homotopy(actualStream(outlet.h_outflow), outlet.h_outflow))
    annotation (Placement(transformation(extent={{66,-10},{86,10}})));
public
  ClaRa.Basics.Interfaces.EyeOut eye if showData annotation (
      Placement(transformation(extent={{90,-70},{110,-50}}),
        iconTransformation(extent={{100,-70},{120,-50}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{45,-61},{47,-59}})));
equation
  eye_int.m_flow = -outlet.m_flow;
  eye_int.T = fluidOut.T-273.15;
  eye_int.s = fluidOut.s/1e3;
  eye_int.p = outlet.p/1e5;
  eye_int.h = fluidOut.h/1e3;

  connect(eye,eye_int)  annotation (Line(
      points={{100,-60},{46,-60}},
      color={255,204,51},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Icon(graphics),
      Diagram(graphics));
end Pump_Base;
