within Exergy.XClaRa.Components.VolumesValvesFittings.Fittings;
model Join_L3_Y "A Y-join with non-ideal mixing"
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

  extends ClaRa.Basics.Icons.Tpipe;

  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");

  outer ClaRa.SimCenter simCenter;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid   medium=simCenter.fluid1
    "Medium in the component"  annotation(Dialog(group="Fundamental Definitions"));

   parameter SI.Volume volume(min=1e-6)=0.1 "System Volume"                               annotation(Dialog(tab="General", group="Geometry"));
  parameter SI.MassFlowRate m_flow_nom= 10 "Nominal mass flow rates at inlet"
                                        annotation(Dialog(tab="General", group="Nominal Values"));
  parameter SI.Pressure p_nom=1e5 "Nominal pressure"  annotation(Dialog(group="Nominal Values"));
//   parameter SI.EnthalpyMassSpecific h_nom=1e5 "Nominal specific enthalpy"  annotation(Dialog(group="Nominal Values"));

  parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"  annotation(Dialog(tab="Initialisation"));
  parameter SI.EnthalpyMassSpecific h_start= 1e5
    "Start value of sytsem specific enthalpy"
                                             annotation(Dialog(tab="Initialisation"));
  parameter SI.Pressure p_start= 1e5 "Start value of sytsem pressure"               annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation" annotation(Dialog(tab="Initialisation"), choicesAllMatching);
    parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "|Summary and Visualisation||True, if expert summary should be applied";
  parameter Boolean showData=true
    "|Summary and Visualisation||True, if a data port containing p,T,h,s,m_flow shall be shown, else false";

  parameter Real y_start=0.5
    "|Initialisation|Start value for relative filling Level";
  parameter SI.EnthalpyMassSpecific h_liq_start=-10 +
      TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(medium,
      p_start) "|Initialisation|Start value of sytsem specific enthalpy";
  parameter SI.EnthalpyMassSpecific h_vap_start=+10 +
      TILMedia.VLEFluidFunctions.dewSpecificEnthalpy_pxi(medium,
      p_start) "|Initialisation|Start value of sytsem specific enthalpy";

  parameter SI.VolumeFraction eps_mix[2]={0.2,0.8}
    "|Mixing Process||Volume fraction V_1/V_tot of min/max mixed outlet";
  replaceable model PressureLoss =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
    annotation (group="Fundamental Definitions",choicesAllMatching=true);
  parameter SI.Time tau_cond=0.03
    "|Mixing Process||Time constant of condensation";
  parameter SI.Time tau_evap=tau_cond
    "|Mixing Process||Time constant of evaporation";
  parameter SI.CoefficientOfHeatTransfer alpha_ph=5000
    "|Mixing Process||HTC of the phase border";
  parameter SI.Area A_phaseBorder=10
    "|Mixing Process||Heat transfer area at phase border";
  parameter Real expHT_phases=0
    "|Mixing Process||Exponent for volume dependency on inter phase HT";

public
  ClaRa.Basics.Interfaces.FluidPortIn inlet2(each Medium=medium)
    "First inlet port" annotation (Placement(transformation(extent={{
            -10,70},{10,90}}), iconTransformation(extent={{-10,90},{
            10,110}})));
public
  ClaRa.Basics.Interfaces.FluidPortIn inlet1(each Medium=medium)
    "First inlet port" annotation (Placement(transformation(extent={{
            -110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.FluidPortOut outlet(Medium=medium) "Outlet port"
                  annotation (Placement(transformation(extent={{90,-10},
            {110,10}})));
public
  ClaRa.Basics.Interfaces.EyeOut eye if showData annotation (
      Placement(transformation(extent={{100,-50},{120,-30}}),
        iconTransformation(extent={{100,-50},{120,-30}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{49,-41},{51,-39}})));
public
  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_L3_TwoZonesNPort
    mixingZone(
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry (
         volume=volume, N_inlet=2),
    medium=medium,
    redeclare model HeatTransfer =
        ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.RealMixed
        (level_rel_start=y_start, eps_mix=eps_mix),
    redeclare model PressureLoss = PressureLoss,
    Tau_cond=tau_cond,
    Tau_evap=tau_evap,
    alpha_ph=alpha_ph,
    A_heat_ph=A_phaseBorder,
    exp_HT_phases=expHT_phases,
    useHomotopy=useHomotopy,
    p_nom=p_nom,
    h_liq_start=h_liq_start,
    h_vap_start=h_vap_start,
    p_start=p_start,
    level_rel_start=y_start,
    initType=initType,
    showExpertSummary=showExpertSummary,
    m_flow_nom=m_flow_nom) annotation (Placement(transformation(
          extent={{-10,-10},{10,10}})));

equation
  connect(mixingZone.inlet[1], inlet1)                 annotation (Line(
      points={{-10,0},{-100,0}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(mixingZone.inlet[2], inlet2)                 annotation (Line(
      points={{-10,0},{-20,0},{-20,58},{0,58},{0,80}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(mixingZone.outlet[1], outlet)                 annotation (Line(
      points={{10,0},{100,0}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));

    eye_int.m_flow=-outlet.m_flow;
    eye_int.T= mixingZone.summary.outlet[1].T-273.15;
    eye_int.s=mixingZone.fluidOut[1].s/1e3;
    eye_int.p=mixingZone.summary.fluid.p[1]/1e5;
    eye_int.h=mixingZone.summary.outlet[1].h/1e3;
  connect(eye_int, eye) annotation (Line(points={{50,-40},{78,-40},{110,-40}}, color={190,190,190}));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=true,  extent={{-100,
            -100},{100,100}})),           Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{100,100}}), graphics));
end Join_L3_Y;
