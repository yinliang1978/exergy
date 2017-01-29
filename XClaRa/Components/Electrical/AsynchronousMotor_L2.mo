within Exergy.XClaRa.Components.Electrical;
model AsynchronousMotor_L2 "A simple asynchronous e-motor"
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
  extends ClaRa.Basics.Icons.Motor;

  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=0,
    powerAux=summary.P_loss) if                                                                                               contributeToCycleSummary;

  import ClaRa.Basics.Units;
  import Modelica.Constants.pi;

  outer ClaRa.SimCenter simCenter;

/////////////////Parameters/////////////////////////////////////
  parameter Boolean activateHeatPort = false
    "True if losses are extracted through heat connector"                                          annotation(Dialog(group="Fundamental Definitions"));
  parameter Units.Power P_nom "Nominal power" annotation(Dialog(group="Fundamental Definitions"));
  parameter Units.Torque tau_bd_nom "Nominal breakdown torque"  annotation(Dialog(group="Fundamental Definitions"));
  parameter Integer N_pp "Number of pole pairs of motor"
                                                        annotation(Dialog(group="Fundamental Definitions"));
  parameter ClaRa.Basics.Units.RPM rpm_nom "nominal speed" annotation(Dialog(group="Fundamental Definitions"));
  parameter Real cosphi(min=0, max=1) "Efficiency factor" annotation(Dialog(group="Fundamental Definitions"));
  parameter Units.Efficiency eta_stator "Lumped constant stator efficiency" annotation(Dialog(group="Fundamental Definitions"));
  parameter Modelica.SIunits.MomentOfInertia J=1500 "Moment of inertia" annotation(Dialog(group="Fundamental Definitions"));

  parameter Units.Frequency f_term_nom = 50 "Nominal excitation frequency" annotation(Dialog(group="Electrics Definitions"));
  parameter Units.ElectricCurrent I_rotor_nom "Rotor nominal current" annotation(Dialog(group="Electrics Definitions"));
  parameter Units.Voltage U_term_nom "Nominal excitation voltage"  annotation(Dialog(group="Electrics Definitions"));

  parameter Boolean useCharLine=false
    "True if characteristic line shall be used (else formula of Kloss)"                                 annotation(Dialog(group="Part Load Definitions"));
  parameter Real charLine_tau_s_[:,2] = [0,2;0.7,1.8;0.95,2.8;1,0]
    "Characteristic line: torque = f(rpm/rpm_nom)"                                                                 annotation(Dialog(group="Part Load Definitions", enable=useCharLine));

  final parameter Real Pi_windings = sqrt(tau_bd_nom*(8*pi^2*f_term_nom^2/N_pp*L_rotor)/3/U_term_nom^2)
    "Ratio rotor windings to stator windings";
  final parameter Units.ElectricResistance R_rotor = P_nom*eta_stator*slip_nom/3/I_rotor_nom^2
    "Electric resistante";
  final parameter Units.Inductance L_rotor =  R_rotor/slip_nom * sqrt(1/cosphi^2-1)/(2*pi*f_term_nom)
    "Inductance of rotor";
  final parameter Real slip_nom = (f_term_nom/N_pp - rpm_nom/60)/(f_term_nom/N_pp)
    "Nominal slip";
  final parameter Units.Torque tau_nom = P_nom/rpm_nom/2/pi*60 "Nominal torque";

  parameter String initOption = "fixed slip" "Init option" annotation(choices(choice="fixed slip", choice="steady state in speed"), Dialog(tab="Initialisation"));
  parameter Real slip_start = (f_term_nom/N_pp - rpm_nom/60)/(f_term_nom/N_pp)
    "Initial slip"                                                                            annotation(Dialog(tab="Initialisation"));

  parameter Boolean showExpertSummary=simCenter.showExpertSummary
    "|Summary and Visualisation||True, if expert summary should be applied";
  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                  annotation(Dialog(tab="Summary and Visualisation"));

  parameter Modelica.Blocks.Types.Smoothness smoothness=Modelica.Blocks.Types.Smoothness.LinearSegments
    "Smoothness of table interpolation"                                                                                                     annotation(Dialog(tab="Expert Settings",enable=useCharLine));

protected
  Units.RPM rpm(start=rpm_nom) "Shaft rotational speed";
  Units.RPM rpm_sync "Synchronous rotational speed";
  Units.Torque tau_bd "Breakdown slip";
  Units.Torque tau_rotor "Air gap electric torque";

  Units.Efficiency eta_gap "Electric efficiency";
  Real slip "Asynchronious motor slip";
  Real slip_bd "Breakdown slip";
  Units.Voltage U_rotor0 "Induced rotor voltage at zero speed";

  Exergy.XClaRa.Components.Utilities.Blocks.ParameterizableTable1D table(
    table=charLine_tau_s_,
    columns={2},
    smoothness=smoothness);

public
model Summary
  extends ClaRa.Basics.Icons.RecordIcon;
  parameter Boolean showExpertSummary
      "True, if expert summary should be applied";
  input Units.Power P_term "Excitation power";
  input Units.Power P_gap "Air gap electric power";
  input Units.Power P_mech "Mechanic power";
  input Units.Power P_loss "Heat loss";
  input Units.RPM rpm "Shaft rotational speed";
  input Units.RPM rpm_sync "Synchronous rotational speed";
  input Real slip "Slip (rpm_sync-rpm)/rpm_sync";
  input Units.Torque tau_rotor "Electric rotor torque";
  input Units.Torque tau_mech "Mechanic shaft torque";
  input Units.Efficiency eta "Total electric efficiency";
  input Units.ElectricVoltage U_term "Terminal volage";
  input Units.ElectricCurrent I_rotor "Rotor current";
end Summary;

/////////////////////////////////////////////////////////
  Modelica.Blocks.Interfaces.RealInput U_term "Terminal voltage" annotation (Placement(transformation(extent={{-140,-20},{-100,20}}), iconTransformation(extent={{-140,-20},{-100,20}})));
  Modelica.Blocks.Interfaces.RealInput f_term "Excitation frequency" annotation (Placement(transformation(
          extent={{-140,20},{-100,60}}),  iconTransformation(extent={{-140,20},{-100,60}})));
  ClaRa.Basics.Interfaces.HeatPort_a    heat(Q_flow= -(summary.P_term-summary.P_mech)) if activateHeatPort
     annotation (Placement(transformation(extent={{10,90},{-10,110}}),
          iconTransformation(extent={{10,90},{-10,110}})));
  Modelica.Mechanics.Rotational.Interfaces.Flange_a shaft "Flange of shaft"
                      annotation (Placement(transformation(extent={{90,-10},{110,10}}), iconTransformation(extent={{90,-10},{110,10}})));
  Summary summary(showExpertSummary=showExpertSummary,
                   P_term = summary.P_gap/eta_stator,
                   P_gap = tau_rotor*der(shaft.phi)/(eta_gap+1e-3),
                   P_mech=-shaft.tau*der(shaft.phi),
                   P_loss=summary.P_term-summary.P_mech,
                   rpm=rpm,
                   rpm_sync=rpm_sync,
                   slip=slip,
                   tau_mech = -shaft.tau,
                   tau_rotor = tau_rotor,
                   eta=eta_stator*eta_gap,
                   U_term=U_term,
                   I_rotor=U_rotor0/sqrt(R_rotor^2/slip^2 + 2*pi*f_term*L_rotor))  annotation (Placement(transformation(
           extent={{-100,-100},{-80,-80}})));

initial equation

  if initOption == "steady state in speed" then
    der(rpm)=0;
  elseif initOption == "fixed slip" then
    slip=slip_start;
  else
    assert(false, "Unknown initOption in e-motor");
  end if;

equation
  //////////////Definitions //////////////////////////
  der(shaft.phi) = (2*pi*rpm/60);
  slip= (rpm_sync-max(0,rpm))/(rpm_sync);
  eta_gap = 1 - slip;
  rpm_sync = f_term*60/N_pp;

///////////////Electric Characteristics /////////////
  U_rotor0 = U_term*Pi_windings;
  tau_bd = 3*U_rotor0^2/(8*pi^2*rpm_sync/60*f_term*L_rotor);
  slip_bd = R_rotor/L_rotor/2/pi/f_term;
  if useCharLine then
    table.u[1] = rpm/rpm_sync;
    tau_rotor = table.y[1]*tau_nom*tau_bd/tau_bd_nom;
  else
    table.u[1] = 1;
    tau_rotor/(tau_bd) = 2/(slip/slip_bd + slip_bd/slip) "Formula of Kloss";
  end if;

///////////////Mechanical Energy Balance ////////////
  J*2*pi/60*rpm*der(rpm) = rpm*(shaft.tau  +tau_rotor) "Energy balance";

  annotation (
    Icon(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                        graphics),
    Diagram(coordinateSystem(preserveAspectRatio=false,extent={{-100,-100},{100,100}}),
                         graphics),
    Documentation(info="<html>
<p>
</html>"));
end AsynchronousMotor_L2;
