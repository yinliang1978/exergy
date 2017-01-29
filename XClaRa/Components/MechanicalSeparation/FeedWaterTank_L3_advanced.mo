within Exergy.XClaRa.Components.MechanicalSeparation;
model FeedWaterTank_L3_advanced
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
  parameter Modelica.SIunits.Length radius_flange=0.05
    "||Geometry|Flange radius";
  parameter SI.Length z_tapping = 0 "||Geometry|position of tapping flange";
  parameter SI.Length z_condensate= 0.1
    "||Geometry|position of condensate flange";
  parameter SI.Length z_aux= 0.1 "||Geometry|position of auxilliary flange";
  parameter SI.Length z_feed = 0 "||Geometry|position of feedwater flange";
  parameter SI.Length z_vent = 0.1 "||Geometry|position of vent flange";

  parameter SI.Time Tau_cond=10 "Time constant of condensation" annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));
  parameter SI.Time Tau_evap=Tau_cond*1000 "Time constant of evaporation" annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));
  parameter Real absorbInflow=1
    "Absorption of incoming mass flow to the zones 1: perfect in the allocated zone, 0: perfect according to steam quality"
                                                                                              annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));
  parameter SI.Area A_phaseBorder=volume.geo.A_hor*100
    "Heat transfer area at phase border"                                                    annotation (Dialog(tab="Phase Separation", group="Heat Transfer Between Phases"));
  parameter SI.CoefficientOfHeatTransfer alpha_ph=500 "HTC of the phase border"
                                                                                annotation (Dialog(tab="Phase Separation", group="Heat Transfer Between Phases"));
  parameter Real expHT_phases=0
    "Exponent for volume dependency on inter phase HT"                             annotation (Dialog(tab="Phase Separation", group="Heat Transfer Between Phases"));
  parameter Boolean equalPressures=true
    "True if pressure in liquid and vapour phase is equal"                                     annotation (Dialog(tab="Phase Separation", group="Mass Transfer Between Phases"));
  parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
    "Smoothness of table interpolation for calculation of filling level"                                annotation(Dialog(tab="Phase Separation", group="Numerical Robustness"));

  parameter SI.EnthalpyMassSpecific h_liq_start=-10 +
      TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(medium,
      volume.p_start)
    "|Initialisation||Start value of liquid specific enthalpy";
  parameter SI.EnthalpyMassSpecific h_vap_start=+10 +
      TILMedia.VLEFluidFunctions.dewSpecificEnthalpy_pxi(medium, volume.p_start)
    "|Initialisation||Start value of vapour specific enthalpy";

  replaceable model PressureLoss =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
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
   input Real level_rel "relative filling level";
 end Outline;

 model Wall
   extends ClaRa.Basics.Icons.RecordIcon;
   input SI.Temperature T_wall[3] "Temperatures";
 end Wall;

public
    outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";
    Exergy.Utilities.ViewObjectNE viewObject(nEnergy={4,0,0,0});

    Exergy.Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{86,88},{106,108}})));

equation
  //  viewObject.h[1].E_flow = inlet.m_flow*noEvent(actualStream(inlet.h_outflow));
  viewObject.h[1].E_flow = summary.tapping.m_flow*summary.tapping.h;
  viewObject.h[1].Ex_flow = summary.tapping.m_flow*(summary.tapping.h-refEnv.T*summary.tapping.s);

    viewObject.h[2].E_flow = summary.condensate.m_flow*summary.condensate.h;
  viewObject.h[2].Ex_flow = summary.condensate.m_flow*(summary.condensate.h-refEnv.T*summary.condensate.s);

    viewObject.h[3].E_flow = summary.aux.m_flow*summary.aux.h;
  viewObject.h[3].Ex_flow = summary.aux.m_flow*(summary.aux.h-refEnv.T*summary.aux.s);

    viewObject.h[4].E_flow = -summary.feedwater.m_flow*summary.feedwater.h;
  viewObject.h[4].Ex_flow = -summary.feedwater.m_flow*(summary.feedwater.h-refEnv.T*summary.feedwater.s);

 // viewObject.E[1].E =liveSteam.d*volume_tot_HP *(liveSteam.h);
 // viewObject.E[1].Ex =liveSteam.d*volume_tot_HP *(liveSteam.h);

 // viewObject.E[1].E =mass_HP *(liveSteam.h)  + mass_IP *(reheatedSteam.h);
 // viewObject.E[1].Ex = mass_HP *(liveSteam.h -  refEnv.T*liveSteam.s)
   //                   + mass_IP *(reheatedSteam.h -  refEnv.T*reheatedSteam.s);
//    viewObject.E[1].E = mass_HP *(liveSteam.h - liveSteam.p/liveSteam.d) + mass_IP *(reheatedSteam.h - reheatedSteam.p/reheatedSteam.d);
  // viewObject.E[1].Ex = mass_HP *(liveSteam.h - liveSteam.p/liveSteam.d-  refEnv.T*liveSteam.s)
 //                     + mass_IP *(reheatedSteam.h - reheatedSteam.p/reheatedSteam.d -  refEnv.T*reheatedSteam.s);

   //viewObject.E[1].Ex = shell.summary.fluid.mass * ( refEnv.T*(shell.summary.fluid.rho));

  connect(viewObject.viewOutput,viewOutput);

public
 model Summary
   extends ClaRa.Basics.Icons.RecordIcon;
   Outline outline;
   Wall wall;
   ClaRa.Basics.Records.FlangeVLE condensate;
   ClaRa.Basics.Records.FlangeVLE tapping;
   ClaRa.Basics.Records.FlangeVLE feedwater;
   ClaRa.Basics.Records.FlangeVLE aux;
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
          volume.fluidIn[1].vleFluidPointer),
      s=TILMedia.VLEFluidObjectFunctions.specificEntropy_phxi(
          heatingSteam.p,
          actualStream(heatingSteam.h_outflow),
          actualStream(heatingSteam.xi_outflow),
          volume.fluidIn[1].vleFluidPointer),
      steamQuality=TILMedia.VLEFluidObjectFunctions.steamMassFraction_phxi(
          heatingSteam.p,
          actualStream(heatingSteam.h_outflow),
          actualStream(heatingSteam.xi_outflow),
          volume.fluidIn[1].vleFluidPointer),
      H_flow=heatingSteam.m_flow*actualStream(heatingSteam.h_outflow),
      rho=TILMedia.VLEFluidObjectFunctions.density_phxi(
          heatingSteam.p,
          actualStream(heatingSteam.h_outflow),
          actualStream(heatingSteam.xi_outflow),
          volume.fluidIn[1].vleFluidPointer)),
    condensate(
      showExpertSummary=showExpertSummary,
      m_flow=condensate.m_flow,
      p=condensate.p,
      h=actualStream(condensate.h_outflow),
      T=TILMedia.VLEFluidObjectFunctions.temperature_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn[2].vleFluidPointer),
      s=TILMedia.VLEFluidObjectFunctions.specificEntropy_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn[2].vleFluidPointer),
      steamQuality=TILMedia.VLEFluidObjectFunctions.steamMassFraction_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn[2].vleFluidPointer),
      H_flow=condensate.m_flow*actualStream(condensate.h_outflow),
      rho=TILMedia.VLEFluidObjectFunctions.density_phxi(
          condensate.p,
          actualStream(condensate.h_outflow),
          actualStream(condensate.xi_outflow),
          volume.fluidIn[2].vleFluidPointer)),
    aux(
      showExpertSummary=showExpertSummary,
      m_flow=aux.m_flow,
      p=aux.p,
      h=actualStream(aux.h_outflow),
      T=TILMedia.VLEFluidObjectFunctions.temperature_phxi(
          aux.p,
          actualStream(aux.h_outflow),
          actualStream(aux.xi_outflow),
          volume.fluidIn[3].vleFluidPointer),
      s=TILMedia.VLEFluidObjectFunctions.specificEntropy_phxi(
          aux.p,
          actualStream(aux.h_outflow),
          actualStream(aux.xi_outflow),
          volume.fluidIn[3].vleFluidPointer),
      steamQuality=TILMedia.VLEFluidObjectFunctions.steamMassFraction_phxi(
          aux.p,
          actualStream(aux.h_outflow),
          actualStream(aux.xi_outflow),
          volume.fluidIn[3].vleFluidPointer),
      H_flow=aux.m_flow*actualStream(aux.h_outflow),
      rho=TILMedia.VLEFluidObjectFunctions.density_phxi(
          aux.p,
          actualStream(aux.h_outflow),
          actualStream(aux.xi_outflow),
          volume.fluidIn[3].vleFluidPointer)),
    feedwater(
      showExpertSummary=showExpertSummary,
      m_flow=-outlet.m_flow,
      p=outlet.p,
      h=actualStream(outlet.h_outflow),
      T=TILMedia.VLEFluidObjectFunctions.temperature_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn[1].vleFluidPointer),
      s=TILMedia.VLEFluidObjectFunctions.specificEntropy_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn[1].vleFluidPointer),
      steamQuality=TILMedia.VLEFluidObjectFunctions.steamMassFraction_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn[1].vleFluidPointer),
      H_flow=-outlet.m_flow*actualStream(outlet.h_outflow),
      rho=TILMedia.VLEFluidObjectFunctions.density_phxi(
          outlet.p,
          actualStream(outlet.h_outflow),
          actualStream(outlet.xi_outflow),
          volume.fluidIn[1].vleFluidPointer))) annotation (Placement(transformation(extent={{-60,-60},{-80,-40}})));
  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLE_L3_TwoZonesNPort volume(
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
    Tau_evap=Tau_evap,
    h_liq_start=h_liq_start,
    h_vap_start=h_vap_start,
    redeclare model Geometry =
        ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry (
        N_heat=1,
        A_heat={Modelica.Constants.pi*diameter*length},
        final A_hor=if orientation == ClaRa.Basics.Choices.GeometryOrientation.vertical
             then Modelica.Constants.pi/4*diameter^2 else diameter*
            length,
        N_outlet=2,
        z_out={z_feed,z_vent},
        shape=if orientation == ClaRa.Basics.Choices.GeometryOrientation.vertical
             then [0,1; 1,1] else [0.0005,0.02981; 0.0245,0.20716;
            0.1245,0.45248; 0.2245,0.58733; 0.3245,0.68065; 0.4245,
            0.74791; 0.5245,0.7954; 0.6245,0.8261; 0.7245,0.84114;
            0.8245,0.84015; 0.9245,0.82031; 1,0.7854],
        height_fill=if orientation == ClaRa.Basics.Choices.GeometryOrientation.vertical
             then length else diameter,
        volume=Modelica.Constants.pi/4*diameter^2*length,
        A_cross=if orientation == ClaRa.Basics.Choices.GeometryOrientation.vertical
             then Modelica.Constants.pi/4*diameter^2 else length*
            diameter,
        final A_front=if orientation == ClaRa.Basics.Choices.GeometryOrientation.vertical
             then Modelica.Constants.pi/4*diameter^2 else length*
            diameter,
        N_inlet=3,
        z_in={z_tapping,z_condensate,z_aux}),
    redeclare model PhaseBorder =
        ClaRa.Basics.ControlVolumes.Fundamentals.SpacialDistribution.RealSeparated
        (
        level_rel_start=level_rel_start,
        radius_flange=radius_flange,
        absorbInflow=absorbInflow,
        smoothness=smoothness),
    A_heat_ph=A_phaseBorder,
    exp_HT_phases=expHT_phases,
    alpha_ph=alpha_ph,
    equalPressures=equalPressures)
    annotation (Placement(transformation(extent={{32,-30},{12,-10}})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.ThickWall_L4 wall(
    sizefunc=+1,
    N_tubes=1,
    initChoice=ClaRa.Basics.Choices.Init.steadyState,
    length=length,
    N_rad=3,
    diameter_i=diameter,
    redeclare replaceable model Material = material,
    diameter_o=diameter + 2*thickness_wall)
             annotation (Placement(transformation(extent={{12,24},{32,44}})));

  ClaRa.Basics.Interfaces.FluidPortOut vent(Medium=medium) annotation (
      Placement(transformation(extent={{-10,88},{10,108}}),
        iconTransformation(extent={{-10,88},{10,108}})));
public
  ClaRa.Basics.Interfaces.EyeOut eye_sat if showData annotation (
      Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-40,110}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=90,
        origin={-40,110})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int1
    annotation (Placement(transformation(extent={{-41,73},{-39,75}})));
public
  ClaRa.Basics.Interfaces.FluidPortIn aux(Medium=medium) "Auxilliary inlet"
                       annotation (Placement(transformation(extent={{
            150,50},{170,70}}), iconTransformation(extent={{150,50},{
            170,70}})));
  Adapters.Scalar2VectorHeatPort scalar2VectorHeatPort(N=2)
    annotation (Placement(transformation(extent={{-10,-10},{10,10}},
        rotation=270,
        origin={22,8})));
  Modelica.Blocks.Interfaces.RealOutput level(value = if outputAbs then summary.outline.level_abs else summary.outline.level_rel) if levelOutput annotation (Placement(transformation(extent={{204,-126},{224,-106}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={242,-110})));
equation
  eye_int1.m_flow=-vent.m_flow;
  eye_int1.T=volume.summary.outlet[2].T-273.15;
  eye_int1.s=volume.fluidOut[2].s/1000;
  eye_int1.h=volume.summary.outlet[2].h/1000;
  eye_int1.p=volume.summary.outlet[2].p/100000;

  eye_int.m_flow=-outlet.m_flow;
  eye_int.T=volume.summary.outlet[1].T-273.15;
  eye_int.s=volume.fluidOut[1].s/1000;
  eye_int.h=volume.summary.outlet[1].h/1000;
  eye_int.p=volume.summary.outlet[1].p/100000;
  connect(volume.inlet[1], heatingSteam) annotation (Line(
      points={{32,-20},{40,-20},{40,60},{-200,60},{-200,80},{-200,80}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.inlet[2], condensate) annotation (Line(
      points={{32,-20},{116,-20},{116,-20},{200,-20},{200,60},{200,60}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.outlet[1], outlet) annotation (Line(
      points={{12,-20},{-260,-20},{-260,-100}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.outlet[2], vent) annotation (Line(
      points={{12,-20},{0,-20},{0,98}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_int1,eye_sat)
                        annotation (Line(
      points={{-40,74},{-40,110}},
      color={190,190,190},
      smooth=Smooth.None));
  connect(aux, volume.inlet[3])           annotation (Line(
      points={{160,60},{44,60},{44,-20},{32,-20}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(scalar2VectorHeatPort.heatVector, volume.heat) annotation (Line(
      points={{22,-2},{22,-10}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(scalar2VectorHeatPort.heatScalar, wall.innerPhase) annotation (Line(
      points={{22,18},{22,24.4},{21.8,24.4}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false,extent={{-300,
            -100},{300,100}}),graphics={
                     Rectangle(extent={{220,-8},{260,-92}}, lineColor={27,36,42},
                               fillColor={153,205,221},
                               fillPattern=FillPattern.Solid,
                               visible=showLevel),
                     Rectangle(extent=DynamicSelect({{220,-50},{260,-92}}, {{220,summary.outline.level_rel*84-92},{260,-92}}),
                               lineColor={27,36,42},
                               fillColor={0,131,169},
                               fillPattern=FillPattern.Solid,
                               visible=showLevel)}),        Diagram(coordinateSystem(
          preserveAspectRatio=true, extent={{-300,-100},{300,100}}),
                                               graphics));
end FeedWaterTank_L3_advanced;
