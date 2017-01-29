within Exergy.XClaRa.Components.MechanicalSeparation;
model FeedWaterTank_L3
  "Feedwater tank : separated volume approach | level-dependent phase separation"
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

extends Exergy.XClaRa.Components.MechanicalSeparation.FeedWaterTank_base;
  parameter ClaRa.Basics.Units.Length thickness_wall=0.005*diameter
    "Thickness of the cylinder wall"                                                                  annotation(Dialog(group="Geometry"));
  replaceable model material = TILMedia.SolidTypes.TILMedia_Steel constrainedby
    TILMedia.SolidTypes.TILMedia_Aluminum "Material of the walls"                                                                              annotation (Dialog(group="Fundamental Definitions"),choicesAllMatching);
  extends ClaRa.Basics.Icons.ComplexityLevel(complexity="L3");
  parameter Modelica.SIunits.Length radius_flange=0.05 "Flange radius" annotation(Dialog(group="Geometry"));
  parameter SI.Time Tau_cond=10 "Time constant of condensation" annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));
  parameter SI.Time tau_evap=Tau_cond*1000 "Time constant of evaporation" annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));
  parameter Real absorbInflow=1
    "Absorption of incoming mass flow to the zones 1: perfect in the allocated zone, 0: perfect according to steam quality"
                                                                                              annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));
  parameter SI.Area A_phaseBorder=volume.geo.A_hor*100
    "Heat transfer area at phase border"                                                    annotation (Dialog(tab="Phase Separation", group="Heat Transfer Between Phases"));
  parameter SI.CoefficientOfHeatTransfer alpha_ph=500 "HTC of the phase border"
                                                                                annotation (Dialog(tab="Phase Separation", group="Heat Transfer Between Phases"));
  //  parameter Real expHT_phases=0 "Exponent for volume dependency on inter phase HT" annotation (Dialog(tab="Phase Separation", group="Heat Transfer Between Phases"));
  parameter Boolean equalPressures=true
    "True if pressure in liquid and vapour phase is equal"                                     annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));

  parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
    "Smoothness of table interpolation for calculation of filling level"                                annotation(Dialog(tab="Phase Separation", group="Numerical Robustness"));
  parameter Modelica.SIunits.Length z_in=1 "Height of inlet ports" annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.Length z_out=1 "Height of outlet ports"  annotation(Dialog(group="Geometry"));
  parameter Modelica.SIunits.SpecificEnthalpy h_liq_start= TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(medium, p_start)
    "|Initialisation||Initial liquid specific enthalpy";
  parameter Modelica.SIunits.SpecificEnthalpy h_vap_start= TILMedia.VLEFluidFunctions.dewSpecificEnthalpy_pxi(medium, p_start)
    "|Initialisation||Initial vapour specific enthalpy";
  replaceable model PressureLoss =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L3
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L3
    "Pressure loss model" annotation(Dialog(group="Fundamental Definitions"), choicesAllMatching);
  replaceable model HeatTransfer =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3
      (                                                                                                    alpha_nom={3000,3000})
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.HeatTransfer_L3
    "Heat transfer to the walls"                                                              annotation (Dialog(group="Fundamental Definitions"),choicesAllMatching=true);

 model Outline
   extends ClaRa.Basics.Icons.RecordIcon;
   input SI.Length level_abs "Absolute filling level";
   input Real level_rel "Relative filling level";
 end Outline;

 model Wall
   extends ClaRa.Basics.Icons.RecordIcon;
   input SI.Temperature T_wall[3] "Temperatures";
 end Wall;

 model Summary
   extends ClaRa.Basics.Icons.RecordIcon;
   Outline outline;
   Wall wall;
   ClaRa.Basics.Records.FlangeVLE condensate;
   ClaRa.Basics.Records.FlangeVLE tapping;
   ClaRa.Basics.Records.FlangeVLE feedwater;
 end Summary;

  Summary summary(
    outline(level_abs=volume.phaseBorder.level_abs, level_rel=volume.phaseBorder.level_rel),
    wall(T_wall=wall.T),
    tapping(
      showExpertSummary=showExpertSummary,
      m_flow=heatingSteam.m_flow,
      p=heatingSteam.p,
      h=actualStream(heatingSteam.h_outflow),
      T=TILMedia.VLEFluidObjectFunctions.temperature_phxi(
          heatingSteam.p,
          actualStream(heatingSteam.h_outflow),
          actualStream(heatingSteam.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      s=TILMedia.VLEFluidObjectFunctions.specificEntropy_phxi(
          heatingSteam.p,
          actualStream(heatingSteam.h_outflow),
          actualStream(heatingSteam.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      steamQuality=TILMedia.VLEFluidObjectFunctions.steamMassFraction_phxi(
          heatingSteam.p,
          actualStream(heatingSteam.h_outflow),
          actualStream(heatingSteam.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      H_flow=heatingSteam.m_flow*actualStream(heatingSteam.h_outflow),
      rho=TILMedia.VLEFluidObjectFunctions.density_phxi(
          heatingSteam.p,
          actualStream(heatingSteam.h_outflow),
          actualStream(heatingSteam.xi_outflow),
          volume.fluidIn.vleFluidPointer)),
    condensate(
      showExpertSummary=showExpertSummary,
      m_flow=condensate.m_flow,
      p=condensate.p,
      h=actualStream(condensate.h_outflow),
      T=TILMedia.VLEFluidObjectFunctions.temperature_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      s=TILMedia.VLEFluidObjectFunctions.specificEntropy_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      steamQuality=TILMedia.VLEFluidObjectFunctions.steamMassFraction_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      H_flow=condensate.m_flow*actualStream(condensate.h_outflow),
      rho=TILMedia.VLEFluidObjectFunctions.density_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn.vleFluidPointer)),
    feedwater(
      showExpertSummary=showExpertSummary,
      m_flow=-outlet.m_flow,
      p=outlet.p,
      h=actualStream(outlet.h_outflow),
      T=TILMedia.VLEFluidObjectFunctions.temperature_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      s=TILMedia.VLEFluidObjectFunctions.specificEntropy_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      steamQuality=TILMedia.VLEFluidObjectFunctions.steamMassFraction_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn.vleFluidPointer),
      H_flow=-outlet.m_flow*actualStream(outlet.h_outflow),
      rho=TILMedia.VLEFluidObjectFunctions.density_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn.vleFluidPointer))) annotation (Placement(transformation(extent={{-60,-60},{-80,-40}})));
  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_3_TwoZones volume(
    medium=medium,
    redeclare model PressureLoss = PressureLoss,
    useHomotopy=useHomotopy,
    m_flow_nom=m_flow_cond_nom + m_flow_heat_nom,
    p_nom=p_nom,
    p_start=p_start,
    initType=initType,
    level_rel_start=level_rel_start,
    Tau_cond=Tau_cond,
    showExpertSummary=showExpertSummary,
    redeclare model HeatTransfer = HeatTransfer,
    h_liq_start=h_liq_start,
    h_vap_start=h_vap_start,
    Tau_evap=tau_evap,
    alpha_ph=alpha_ph,
    A_heat_ph=A_phaseBorder,
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.RealSeparated
        (
        level_rel_start=level_rel_start,
        smoothness=smoothness,
        radius_flange=radius_flange,
        absorbInflow=absorbInflow),
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowCylinder (
        orientation=orientation,
        diameter=diameter,
        length=length,
        z_in={z_in},
        z_out={z_out}),
    equalPressures=equalPressures)
    annotation (Placement(transformation(extent={{12,-30},{-8,-10}})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.ThickWall_L4 wall(
    sizefunc=+1,
    N_tubes=1,
    initChoice=ClaRa.Basics.Choices.Init.steadyState,
    length=length,
    N_rad=3,
    diameter_i=diameter,
    redeclare replaceable model Material = material,
    diameter_o=diameter + 2*thickness_wall)
             annotation (Placement(transformation(extent={{-8,6},{12,26}})));

  Modelica.Blocks.Interfaces.RealOutput level(value = if outputAbs then summary.outline.level_abs else summary.outline.level_rel) if levelOutput annotation (Placement(transformation(extent={{204,-126},{224,-106}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={240,-110})));
equation
  eye_int.m_flow=-outlet.m_flow;
  eye_int.T=volume.summary.outlet.T-273.15;
  eye_int.s=volume.fluidOut.s/1000;
  eye_int.h=volume.summary.outlet.h/1000;
  eye_int.p=volume.summary.outlet.p/100000;

  connect(volume.inlet, condensate) annotation (Line(
      points={{12,-20},{106,-20},{106,60},{200,60}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(heatingSteam, volume.inlet) annotation (Line(
      points={{-200,80},{-200,52},{28,52},{28,-20},{12,-20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.outlet, outlet) annotation (Line(
      points={{-8,-20},{-134,-20},{-134,-100},{-260,-100}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(wall.innerPhase, volume.heat) annotation (Line(
      points={{1.8,6.4},{2,6.4},{2,-10}},
      color={191,0,0},
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=true, extent={{-300,
            -100},{300,100}}),
                   graphics={
                     Rectangle(extent={{220,-8},{260,-92}}, lineColor={27,36,42},
                               fillColor={153,205,221},
                               fillPattern=FillPattern.Solid,
                               visible=showLevel),
                     Rectangle(extent=DynamicSelect({{220,-50},{260,-92}}, {{220,summary.outline.level_rel*84-92},{260,-92}}),
                               lineColor={27,36,42},
                               fillColor={0,131,169},
                               fillPattern=FillPattern.Solid,
                               visible=showLevel)}),          Diagram(coordinateSystem(
          preserveAspectRatio=true, extent={{-220,-120},{120,100}})));
end FeedWaterTank_L3;
