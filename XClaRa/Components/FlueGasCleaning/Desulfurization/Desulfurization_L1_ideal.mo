within Exergy.XClaRa.Components.FlueGasCleaning.Desulfurization;
model Desulfurization_L1_ideal
  "Model for an idealised desulfurization with chalk washing"
//___________________________________________________________________________//
// Package of the ClaRa library, version: 1.1.2                              //
// Models of the ClaRa library are tested under DYMOLA v2016 FD01.           //
// It is planned to support alternative Simulators like SimulationX in the   //
// future                                                                    //
//___________________________________________________________________________//
// Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
// Copyright © 2013-2016, DYNCAP/DYNSTART research team.                     //
//___________________________________________________________________________//
// This Modelica package is free software and the use is completely at your  //
// own risk; it can be redistributed and/or modified under the terms of the  //
// Modelica License 2. For license conditions (including the disclaimer of   //
// warranty) see Modelica.UsersGuide.ModelicaLicense2 or visit               //
// http://www.modelica.org/licenses/ModelicaLicense2                         //
//___________________________________________________________________________//
// DYNCAP and DYNSTART are research projects supported by the German Federal //
// Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
// The research team consists of the following project partners:             //
// Institute of Energy Systems (Hamburg University of Technology),           //
// Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
// TLK-Thermo GmbH (Braunschweig, Germany),                                  //
// XRG Simulation GmbH (Hamburg, Germany).                                   //
//___________________________________________________________________________//
  outer ClaRa.SimCenter simCenter;
  extends ClaRa.Basics.Icons.Separator;

  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=-P_el,
    powerAux=0) if contributeToCycleSummary;

  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=simCenter.flueGasModel)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortOut outlet(Medium=simCenter.flueGasModel)
    annotation (Placement(transformation(extent={{90,-10},{110,10}}),
        iconTransformation(extent={{90,-10},{110,10}})));

//## S U M M A R Y   D E F I N I T I O N ###################################################################
   model Outline
    extends ClaRa.Basics.Icons.RecordIcon;
    input Modelica.SIunits.Volume volume "System volume"
      annotation (Dialog(show));
    input Modelica.SIunits.Mass m "System mass" annotation (Dialog(show));
    input Modelica.SIunits.Enthalpy H "System enthalpy"
      annotation (Dialog(show));
    input Modelica.SIunits.Pressure p "System pressure"
      annotation (Dialog(show));
    input Modelica.SIunits.Pressure Delta_p "Pressure loss"      annotation (Dialog(show));
    input Modelica.SIunits.SpecificEnthalpy h "System specific enthalpy"        annotation(Dialog(show));

    input Modelica.SIunits.Temperature T "System temperature"
      annotation (Dialog(show));
    input Modelica.SIunits.Power P_el "Electric power consumption"
      annotation (Dialog(show));
    input Modelica.SIunits.MassFlowRate m_flow_SOx "Separated SOx flow rate"
      annotation (Dialog(show));
    input Real SOx_separationRate "NOx separation rate"
      annotation (Dialog(show));
    input Modelica.SIunits.MassFlowRate m_flow_CaCO3 "Required CaCO3 flow rate"
      annotation (Dialog(show));
    input Modelica.SIunits.MassFlowRate m_flow_CaSO4_H2O
      "Outlet CaSO4_H2O flow rate"
      annotation (Dialog(show));
    input Modelica.SIunits.MassFlowRate m_flow_H2O "Required H2O flow rate"
      annotation (Dialog(show));

   end Outline;

 model Summary
     extends ClaRa.Basics.Icons.RecordIcon;
     Outline outline;
     ClaRa.Basics.Records.FlangeGas  inlet;
     ClaRa.Basics.Records.FlangeGas  outlet;
 end Summary;

//_____________defintion of medium used in cell__________________________________________________________
  inner parameter TILMedia.GasTypes.BaseGas      medium = simCenter.flueGasModel
    "Medium to be used in tubes" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  parameter Real SOx_separationRate = 0.95 "Sulphur separation rate" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));
  parameter ClaRa.Basics.Units.Temperature T_in_H2O = 313.15
    "Temperature of water inlet"                                                          annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));
  parameter Real specificPowerConsumption(unit="J/m3") = 9000
    "Specific power consumption per standard m^3"                                                           annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  replaceable model PressureLoss =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L2
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching=true);

  replaceable model Geometry =
      ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.HollowCylinder(diameter=4,length=10,z_in={0},z_out={10},orientation = ClaRa.Basics.Choices.GeometryOrientation.vertical,flowOrientation = ClaRa.Basics.Choices.GeometryOrientation.vertical)
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Geometry"), choicesAllMatching=true);

  inner parameter Modelica.SIunits.MassFlowRate m_flow_nom= 200
    "Nominal mass flow rates at inlet"                                                             annotation(Dialog(tab="General", group="Nominal Values"));
  inner parameter ClaRa.Basics.Units.Pressure p_nom=1e5 "Nominal pressure"                    annotation(Dialog(group="Nominal Values"));
  inner parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_nom=1e5
    "Nominal specific enthalpy"                                                                annotation(Dialog(group="Nominal Values"));

  inner parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation"                                                                                   annotation(Dialog(tab="Initialisation", choicesAllMatching));
  parameter ClaRa.Basics.Units.Temperature T_start= 273.15 + 100.0
    "Start value of system temperature"                                                                  annotation(Dialog(tab="Initialisation"));

  parameter ClaRa.Basics.Units.Pressure p_start= 1.013e5
    "Start value of sytsem pressure"                                                      annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.MassFraction xi_start[medium.nc-1]={0.01,0,0.25,0,0.7,0,0,0.04,0}
    "Start value of system mass fraction"                                                                                              annotation(Dialog(tab="Initialisation"));
  inner parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"                                                         annotation(Dialog(tab="Initialisation"));

  parameter Boolean allow_reverseFlow = true
    "True if simulation shall stop at reverse flow conditions"                                          annotation(Dialog(tab="Expert Settings", group="General"));
  parameter Boolean use_dynamicMassbalance = true
    "True if species balance shall be dynamic"                                                annotation(Dialog(tab="Expert Settings", group="General"));

  parameter Boolean useStabilisedMassFlow=false
    "True if the outlet mass flow shall be low-pass filtered"                                             annotation(Dialog(tab="Expert Settings", group="Numerical Robustness"));
  parameter SI.Time Tau= 0.001 "Time Constant of Stabilisation" annotation(Dialog(tab="Expert Settings", group = "Numerical Robustness", enable=useStabilisedMassFlow));

  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));

ClaRa.Basics.Units.Power P_el "Electric power consumption";
ClaRa.Basics.Units.VolumeFlowRate V_flow_std "Standardized volume flow rate";

    TILMedia.GasObjectFunctions.GasPointer GasPointer=
      TILMedia.GasObjectFunctions.GasPointer(
      medium.concatGasName,
      8,
      medium.xi_default,
      medium.nc_propertyCalculation,
      medium.nc,
      0,
      0) "Pointer to external medium memory";

  Fundamentals.Desulfurisation_controlVolume_ideal deSO_controlVolume(
    useStabilisedMassFlow=useStabilisedMassFlow,
    Tau=Tau,
    SOx_separationRate=SOx_separationRate,
    T_in_H2O=T_in_H2O,
    xi_start=xi_start) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={-34,0})));

  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell(
    redeclare model Geometry = Geometry,
    redeclare model PressureLoss = PressureLoss,
    T_start=T_start,
    p_start=p_start,
    xi_start=xi_start,
    initType=initType,
    m_flow_nom=m_flow_nom,
    p_nom=p_nom,
    h_nom=h_nom,
    useHomotopy=useHomotopy,
    allow_reverseFlow=allow_reverseFlow,
    use_dynamicMassbalance=use_dynamicMassbalance,
    redeclare model HeatTransfer =
        Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.IdealHeatTransfer_L2)
    annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={26,0})));

inner Summary summary(outline(
    volume=flueGasCell.summary.outline.volume_tot,
    m=flueGasCell.summary.outline.mass,
    H=flueGasCell.summary.outline.H,
    h=flueGasCell.summary.outline.h,
    T=flueGasCell.summary.outline.T,
    p=flueGasCell.summary.outline.p,
    Delta_p=flueGasCell.pressureLoss.Delta_p,
    P_el = P_el,
    SOx_separationRate = SOx_separationRate,
    m_flow_SOx = deSO_controlVolume.m_flow_SOx_sep,
    m_flow_CaCO3 = deSO_controlVolume.m_flow_CaCO3_req,
    m_flow_CaSO4_H2O = deSO_controlVolume.m_flow_CaSO4_H2O_out,
    m_flow_H2O = deSO_controlVolume.m_flow_H2O_req),
    inlet(m_flow = inlet.m_flow,
          T = inStream(inlet.T_outflow),
          p = inlet.p,
          h = deSO_controlVolume.flueGasInlet.h,
          xi = inStream(inlet.xi_outflow),
          H_flow = inlet.m_flow*deSO_controlVolume.flueGasInlet.h),
    outlet(m_flow = -outlet.m_flow,
          T = outlet.T_outflow,
          p = inlet.p,
          h = flueGasCell.flueGasOutlet.h,
          xi = outlet.xi_outflow,
          H_flow = -outlet.m_flow*flueGasCell.flueGasOutlet.h)) annotation (Placement(transformation(extent={{-100,-100},{-80,-80}})));

public
  ClaRa.Basics.Interfaces.EyeOut eyeOut annotation (Placement(
        transformation(extent={{80,-78},{120,-42}}),
        iconTransformation(extent={{90,-50},{110,-30}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{32,-68},{48,-52}}),
        iconTransformation(extent={{90,-84},{84,-78}})));
equation
  V_flow_std = inlet.m_flow / TILMedia.GasObjectFunctions.density_pTxi(1.01325e5,273.15,inStream(inlet.xi_outflow),GasPointer);
  P_el = specificPowerConsumption * V_flow_std;

  //______________Eye port variable definition________________________
  eye_int.m_flow = -outlet.m_flow;
  eye_int.T = flueGasCell.bulk.T-273.15;
  eye_int.s = flueGasCell.bulk.s/1e3;
  eye_int.p = flueGasCell.bulk.p/1e5;
  eye_int.h = flueGasCell.bulk.h/1e3;

  connect(inlet, deSO_controlVolume.inlet) annotation (Line(
      points={{-100,0},{-44,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(deSO_controlVolume.outlet, flueGasCell.inlet) annotation (Line(
      points={{-24,0},{16,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasCell.outlet, outlet) annotation (Line(
      points={{36,0},{100,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_int,eyeOut)  annotation (Line(
      points={{40,-60},{100,-60}},
      color={190,190,190},
      smooth=Smooth.None));
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,
            -100},{100,100}}),
                         graphics={                  Text(
          extent={{-84,72},{-12,28}},
          lineColor={27,36,42},
          textString="SOx")}),   Diagram(coordinateSystem(preserveAspectRatio=false,
                   extent={{-100,-100},{100,100}})),
    Documentation(info="<html>
<p><b>Information</b></p>
<p><b>Model description: </b>An ideal desulfurization model</p>
<p><b>Contact:</b> Andre Th&uuml;ring, Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>This model uses TILMedia</li>
<li>Calculates the separated SO2 with a given separation rate</li>
<li>Power consumption is calculated with a given specific power consumption </li>
<li>Stationary mass and energy balance inside the ideal control volume</li>
<li>Dwelltime is regarded by upstream flue gas cell</li>
<li>The model is adiabatic</li>
<li>The flue gas leaves the control volume saturated with water</li>
<li>Outlet temperature is calculated with evaporation heat of spray water needed to saturate the flue gas</li>
</ul></p>
</html>"));
end Desulfurization_L1_ideal;
