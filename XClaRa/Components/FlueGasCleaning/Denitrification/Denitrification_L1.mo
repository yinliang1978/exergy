within Exergy.XClaRa.Components.FlueGasCleaning.Denitrification;
model Denitrification_L1
  "Model for a simple ammonia denitrification with fixed separation ratio"
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

extends ClaRa.Basics.Icons.Separator;
  outer ClaRa.SimCenter simCenter;

//## S U M M A R Y   D E F I N I T I O N ###################################################################
 model Outline
    extends ClaRa.Basics.Icons.RecordIcon;
    input Modelica.SIunits.Volume volume "System volume"
      annotation (Dialog(show));
    input Modelica.SIunits.Mass mass "System mass" annotation (Dialog(show));
    input Modelica.SIunits.Enthalpy H "System enthalpy"
      annotation (Dialog(show));
    input Modelica.SIunits.Pressure p "System pressure"
      annotation (Dialog(show));
    input Modelica.SIunits.Pressure Delta_p "Pressure loss"      annotation (Dialog(show));
    input Modelica.SIunits.SpecificEnthalpy h "System specific enthalpy"        annotation(Dialog(show));

    input Modelica.SIunits.Temperature T "System temperature"
      annotation (Dialog(show));
    input Modelica.SIunits.MassFlowRate m_flow_NH3 "Requirered NH3 flow rate"
      annotation (Dialog(show));
    input Modelica.SIunits.MassFlowRate m_flow_O2 "Requirered O2 flow rate"
      annotation (Dialog(show));
    input Modelica.SIunits.HeatFlowRate reactionHeat
      "Reaction heat of deNOx catalysis" annotation (Dialog(show));
    input Real NOx_separationRate "NOx separation rate"
      annotation (Dialog(show));
 end Outline;

 model Summary
     extends ClaRa.Basics.Icons.RecordIcon;
     Outline outline;
     ClaRa.Basics.Records.FlangeGas  inlet;
     ClaRa.Basics.Records.FlangeGas  outlet;
 end Summary;

//## P A R A M E T E R S #######################################################################################
//_____________defintion of medium used in cell__________________________________________________________
  inner parameter TILMedia.GasTypes.BaseGas      medium = simCenter.flueGasModel
    "Medium to be used in tubes"                                                                              annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  parameter Real separationRate(max = 0.99995) = 0.9995 "Separation rate" annotation (Dialog(group="Fundamental Definitions"));

  inner parameter ClaRa.Basics.Units.MassFlowRate m_flow_nom= 200
    "Nominal mass flow rates at inlet"                                                                annotation(Dialog(tab="General", group="Nominal Values"));
  inner parameter ClaRa.Basics.Units.Pressure p_nom=1e5 "Nominal pressure"                    annotation(Dialog(group="Nominal Values"));
  inner parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_nom=1e5
    "Nominal specific enthalpy"                                                                annotation(Dialog(group="Nominal Values"));

  inner parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation"                                                                                   annotation(Dialog(tab="Initialisation", choicesAllMatching));
  parameter ClaRa.Basics.Units.Temperature T_start= 273.15 + 100.0
    "Start value of system temperature"                                                                annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.Pressure p_start= 1.013e5
    "Start value of sytsem pressure"                                                      annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.MassFraction xi_start[medium.nc-1]={0.01,0,0.25,0,0.7,0,0,0.04,0}
    "Start value of sytsem mass fraction"                                                                                              annotation(Dialog(tab="Initialisation"));
  inner parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"                                                         annotation(Dialog(tab="Initialisation"));

  parameter Boolean allow_reverseFlow = true
    "True if flow reversal shall be supported"                                          annotation(Evaluate=true, Dialog(tab="Expert Settings"));
  parameter Boolean use_dynamicMassbalance = true
    "True if a dynamic mass balance shall be applied"                                               annotation(Evaluate=true, Dialog(tab="Expert Settings"));

  parameter Boolean showData=true
    "True, if a data port containing p,T,h,s,m_flow shall be shown, else false"
                                                                                            annotation (Dialog(tab="Summary and Visualisation"));
//## V A R I A B L E   P A R T##################################################################################

//____Connectors________________________________________________________________________________________________
  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium, m_flow(min=if
          allow_reverseFlow then -Modelica.Constants.inf else 1e-5))
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortOut outlet(Medium=medium, m_flow(max=if
          allow_reverseFlow then Modelica.Constants.inf else -1e-5))
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat annotation (Placement(transformation(extent={{-10,90},
            {10,110}}),        iconTransformation(extent={{-62,86},{-42,106}})));

//____replaceable models for heat transfer, pressure loss and geometry_____________________________________________________________________________________
replaceable model Geometry =
      ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Geometry"), choicesAllMatching=true);
 replaceable model HeatTransfer =
      ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.Adiabat_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.HeatTransport.Generic_HT.HeatTransfer_L2
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching=true);
  replaceable model PressureLoss =
      ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.NoFriction_L2
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.PressureLoss.Generic_PL.PressureLoss_L2
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Fundamental Definitions"), choicesAllMatching=true);

  Exergy.XClaRa.Components.FlueGasCleaning.Denitrification.Fundamentals.Denitrification_controlVolume
    deNOx_controlVolume(separationRate=separationRate,
      T_NH3_O2_mixture=473.15) annotation (Placement(transformation(
          extent={{-64,-20},{-24,20}})));
  ClaRa.Basics.ControlVolumes.GasVolumes.VolumeGas_L2 flueGasCell(
    redeclare model Geometry = Geometry,
    redeclare model HeatTransfer = HeatTransfer (heatSurfaceAlloc=1),
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
    use_dynamicMassbalance=use_dynamicMassbalance) annotation (Placement(transformation(extent={{12,-22},{56,22}})));

inner Summary summary(outline(
    volume=flueGasCell.summary.outline.volume_tot,
    mass=flueGasCell.summary.outline.mass,
    H=flueGasCell.summary.outline.H,
    h=flueGasCell.summary.outline.h,
    T=flueGasCell.summary.outline.T,
    p=flueGasCell.summary.outline.p,
    Delta_p=flueGasCell.pressureLoss.Delta_p,
    NOx_separationRate = separationRate,
    m_flow_NH3 = deNOx_controlVolume.n_flow_NH3_req*deNOx_controlVolume.NH3_O2_in.M_i[
        9],
    m_flow_O2 = deNOx_controlVolume.n_flow_O2_req*deNOx_controlVolume.NH3_O2_in.M_i[
        6],
    reactionHeat=deNOx_controlVolume.Qdot),
    inlet(m_flow = inlet.m_flow,
          T = inStream(inlet.T_outflow),
          p = inlet.p,
          h = deNOx_controlVolume.flueGasInlet.h,
          xi = inStream(inlet.xi_outflow),
          H_flow = deNOx_controlVolume.flueGasInlet.h*inlet.m_flow),
    outlet(m_flow = -outlet.m_flow,
          T = outlet.T_outflow,
          p = inlet.p,
          h = flueGasCell.flueGasOutlet.h,
          xi = outlet.xi_outflow,
          H_flow = -flueGasCell.flueGasOutlet.h*outlet.m_flow)) annotation (Placement(transformation(extent={{-100,
            -114},{-80,-94}})));

public
  ClaRa.Basics.Interfaces.EyeOut eyeOut if showData annotation (
      Placement(transformation(extent={{80,-78},{120,-42}}),
        iconTransformation(extent={{90,-50},{110,-30}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{48,-68},{32,-52}}),
        iconTransformation(extent={{90,-84},{84,-78}})));

equation
  //______________Eye port variable definition________________________
  eye_int.m_flow = -outlet.m_flow;
  eye_int.T = flueGasCell.bulk.T-273.15;
  eye_int.s = flueGasCell.bulk.s/1e3;
  eye_int.p = flueGasCell.bulk.p/1e5;
  eye_int.h = flueGasCell.bulk.h/1e3;

  connect(inlet, deNOx_controlVolume.inlet) annotation (Line(
      points={{-100,0},{-64,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(deNOx_controlVolume.outlet, flueGasCell.inlet) annotation (Line(
      points={{-24,0},{12,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(flueGasCell.outlet, outlet) annotation (Line(
      points={{56,0},{100,0}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  connect(heat, flueGasCell.heat) annotation (Line(
      points={{0,100},{0,46},{34,46},{34,22},{34,22}},
      color={167,25,48},
      thickness=0.5,
      smooth=Smooth.None));
  connect(eye_int,eyeOut)  annotation (Line(
      points={{40,-60},{100,-60}},
      color={190,190,190},
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,
            -100},{100,100}}),
                      graphics), Icon(coordinateSystem(preserveAspectRatio=false,
          extent={{-100,-100},{100,100}}), graphics={Text(
          extent={{-84,72},{-12,28}},
          lineColor={27,36,42},
          textString="NOx")}),
    Documentation(info="<html>
<p><b>Model description: </b>A simple DeNOx filter model</p>
<p><b>Contact:</b> Andre Th&uuml;ring, Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>This model uses TILMedia</li>
<li>Calculates the amount on NH3 necessary to achieve the given separation rate</li>
<li>Stationary mass and energy balance</li>
</ul></p>
</html>"));
end Denitrification_L1;
