within Exergy.XClaRa.Components.MechanicalSeparation;
model BalanceTank_L3 "A balance tank with a vent"
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
 extends ClaRa.Basics.Icons.BalanceTank;
  outer ClaRa.SimCenter simCenter;

  parameter TILMedia.VLEFluidTypes.BaseVLEFluid liquidMedium=simCenter.fluid1
    "Liquid medium in the component"                                                                           annotation (choicesAllMatching=true,Dialog(group="Fundamental Definitions"));
  parameter TILMedia.GasTypes.BaseGas gasMedium=simCenter.flueGasModel
    "Gas medium in the component"                                                                    annotation (choicesAllMatching=true,Dialog(group="Fundamental Definitions"));
  replaceable model Material = TILMedia.SolidTypes.TILMedia_Aluminum constrainedby
    TILMedia.SolidTypes.BaseSolid "Solid material of the tank"
    annotation (choicesAllMatching=true,Dialog(group="Fundamental Definitions"));
  replaceable model HeatTransfer =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Constant_L3
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.HeatTransfer_L3
    "Heat transfer model"
    annotation (choicesAllMatching=true,Dialog(group="Fundamental Definitions"));
  replaceable model PressureLoss =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.LinearParallelZones_L3
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L3
    "|Fundamental Definitions|Pressure loss model"
    annotation (choicesAllMatching=true,Dialog(group="Fundamental Definitions"));

  parameter ClaRa.Basics.Units.Length diameter_i=2
    "|Geometry|Inner diameter of the tank";
  parameter ClaRa.Basics.Units.Length s_wall=2
    "|Geometry|Wall thickness of the tank";
  parameter ClaRa.Basics.Units.Length height=2 "|Geometry|Height of the tank";
  parameter ClaRa.Basics.Units.Temperature T_start[3]=ones(3)*293.15
    "|Initialisation|Wall|Start values of wall temperature";
  parameter ClaRa.Basics.Choices.Init initWall=ClaRa.Basics.Choices.Init.noInit
    "|Initialisation|Wall|Wall init option"
                                           annotation (choicesAllMatching);

  parameter ClaRa.Basics.Units.Length z_in[3]=ones(3)*height
    "|Geometry|Height of liquid inlet ports";
  parameter ClaRa.Basics.Units.Length z_out[1]={0.1}
    "|Geometry|Height of liquid outlet ports";

  parameter ClaRa.Basics.Units.CoefficientOfHeatTransfer alpha_ph=500
    "|Expert Settings|Phase Border|HTC of the phase border";
  parameter ClaRa.Basics.Units.Area A_phaseBorder=volume.geo.A_hor*100
    "|Expert Settings|Phase Border|Heat transfer area at phase border";

  parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_liq_start=-10 +
      TILMedia.VLEFluidFunctions.bubbleSpecificEnthalpy_pxi(liquidMedium,
      p_start) "|Initialisation|Fluids|Start value ofliquid specific enthalpy";
  parameter ClaRa.Basics.Units.Temperature T_gas_start=293.15
    "|Initialisation|Fluids|Start value of gas zone's temperature";
  parameter ClaRa.Basics.Units.Pressure p_start=1e5
    "|Initialisation|Fluids|Start value of sytsem pressure";
  parameter ClaRa.Basics.Units.MassFraction xi_start[gasMedium.nc - 1]=zeros(gasMedium.nc
       - 1) "|Initialisation|Fluids|Initial gas mass fraction";
  parameter Real relLevel_start=0.5
    "|Initialisation|Fluids|Initial value for relative level";
  parameter String initFluid="No init, use start values as guess"
    "|Initialisation|Fluids|Type of initialisation" annotation (choices(choice = "No init, use start values as guess", choice="Steady state in p, h_liq, T_gas",
            choice = "Steady state in p", choice="steady State in h_liq and T_gas", choice = "Fixed value for filling level",
             choice = "Fixed values for filling level, p, h_liq, T_gas"));

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "True, if expert summary should be applied"                                                               annotation(Dialog(enable = levelOutput, tab="Summary and Visualisation"));
  parameter Boolean levelOutput = false
    "True, if Real level connector shall be addded"                                      annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean outputAbs = false "True, if absolute level is at output"  annotation(Dialog(enable = levelOutput, tab="Summary and Visualisation"));

  ClaRa.Basics.Interfaces.FluidPortOut outlet(Medium=liquidMedium)
    "Outlet port" annotation (Placement(transformation(extent={{-106,-68},
            {-86,-48}}), iconTransformation(extent={{-106,-68},{-86,-48}})));
  ClaRa.Basics.Interfaces.FluidPortIn inlet3(Medium=liquidMedium) "Inlet port"
                 annotation (Placement(transformation(extent={{170,190},
            {190,210}}), iconTransformation(extent={{170,190},{190,210}})));
  ClaRa.Basics.Interfaces.FluidPortIn inlet1(Medium=liquidMedium) "Inlet port"
                 annotation (Placement(transformation(extent={{90,190},
            {110,210}}), iconTransformation(extent={{90,190},{110,210}})));
  ClaRa.Basics.Interfaces.FluidPortIn inlet2(Medium=liquidMedium) "Inlet port"
                 annotation (Placement(transformation(extent={{132,190},
            {152,210}}), iconTransformation(extent={{132,190},{152,210}})));

  ClaRa.Basics.ControlVolumes.SolidVolumes.ThickWall_L4 wall(
    N_rad=3,
    sizefunc=+1,
    diameter_i=diameter_i,
    length=height,
    T_start=T_start,
    diameter_o=diameter_i + 2*s_wall,
    redeclare model Material = Material,
    initChoice=initWall)
    annotation (Placement(transformation(extent={{32,-36},{52,-16}})));

  ClaRa.Basics.ControlVolumes.FluidVolumes.VolumeVLEGas_L3 volume(
    medium=liquidMedium,
    gasType=gasMedium,
    redeclare model HeatTransfer = HeatTransfer,
    redeclare model PressureLoss = PressureLoss,
    alpha_ph=alpha_ph,
    A_heat_ph=A_phaseBorder,
    h_liq_start=h_liq_start,
    T_gas_start=T_gas_start,
    p_start=p_start,
    xi_start=xi_start,
    showExpertSummary=showExpertSummary,
    level_rel_start=relLevel_start,
    initType=initFluid,
    redeclare model Geometry =
        Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry (
        volume=height*diameter_i^2*Modelica.Constants.pi/4,
        N_heat=1,
        A_heat={Modelica.Constants.pi*diameter_i*height},
        A_cross=diameter_i^2*Modelica.Constants.pi/4,
        A_front=diameter_i^2*Modelica.Constants.pi/4,
        A_hor=diameter_i^2*Modelica.Constants.pi/4,
        N_inlet=3,
        N_outlet=1,
        height_fill=height,
        shape=[0,1; 1,1],
        z_out=z_out,
        z_in=z_in))
    annotation (Placement(transformation(extent={{52,-68},{32,-48}})));

  ClaRa.Basics.Interfaces.GasPortIn vent1(Medium=gasMedium) annotation (
     Placement(transformation(extent={{40,190},{60,210}}),
        iconTransformation(extent={{40,190},{60,210}})));

  Modelica.Blocks.Interfaces.RealOutput level(value = if outputAbs then volume.level_abs else volume.level_rel) if levelOutput annotation (Placement(transformation(extent={{204,-126},{224,-106}}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={160,-110})));
equation
  connect(volume.outlet[1], outlet) annotation (Line(
      points={{32,-58},{-96,-58}},
      color={0,131,169},
      pattern=LinePattern.Solid,
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.inlet[1], inlet1) annotation (Line(
      points={{52,-58},{100,-58},{100,200}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.inlet[2], inlet2) annotation (Line(
      points={{52,-58},{142,-58},{142,200}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.inlet[3], inlet3) annotation (Line(
      points={{52,-58},{180,-58},{180,200}},
      color={0,131,169},
      thickness=0.5,
      smooth=Smooth.None));
  connect(volume.vent, vent1) annotation (Line(
      points={{52,-54.2},{52,200},{50,200}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(wall.innerPhase, volume.heat) annotation (Line(
      points={{41.8,-35.6},{42,-35.6},{42,-48.2}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{200,200}}), graphics), Icon(coordinateSystem(
          preserveAspectRatio=false, extent={{-100,-100},{200,200}})));
end BalanceTank_L3;
