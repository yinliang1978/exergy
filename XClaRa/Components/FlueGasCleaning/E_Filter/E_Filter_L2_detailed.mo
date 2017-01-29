within Exergy.XClaRa.Components.FlueGasCleaning.E_Filter;
model E_Filter_L2_detailed
  "Model for an electrical dust filter based on the Deutsch-Equation for the separation rate, migration speed calculated in accodance with  C. Riehle, Basic and Theoretical Operation of ESPs, (1997)"
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

 extends ClaRa.Basics.Icons.E_Filter;
  outer ClaRa.SimCenter simCenter;
  ClaRa.Basics.Interfaces.Connected2SimCenter connected2SimCenter(
    powerIn=0,
    powerOut=-powerConsumption,
    powerAux=0) if contributeToCycleSummary;

//## S U M M A R Y   D E F I N I T I O N ###################################################################
 model Outline
  //  parameter Boolean showExpertSummary annotation(Dialog(hide));
    input Modelica.SIunits.Volume V "System volume"      annotation (Dialog(show));
    input Modelica.SIunits.Mass m "System mass" annotation (Dialog(show));
    input Modelica.SIunits.Enthalpy H "System enthalpy"      annotation (Dialog(show));
    input Modelica.SIunits.SpecificEnthalpy h "System specific enthalpy"        annotation(Dialog(show));
    input Modelica.SIunits.Pressure p "System pressure"      annotation (Dialog(show));
    input Modelica.SIunits.Pressure Delta_p "Pressure loss"      annotation (Dialog(show));
    input Modelica.SIunits.Length lambda "Mean free path of particles";
    input Real Cu "Cunningham slip correction factor"     annotation (Dialog(show));
    input Modelica.SIunits.Velocity w_m
      "Migration speed of dust particles in the E-field"                                  annotation (Dialog(show));
    input Modelica.SIunits.ElectricFieldStrength E_applied "E-Field in Filter" annotation (Dialog(show));
    input Modelica.SIunits.ElectricCharge Q_sat
      "Saturation charge of particles"                                           annotation (Dialog(show));
    input Real separationRate "Separation rate of E-Filter";
    input Modelica.SIunits.MassFlowRate m_flow_dust_out
      "Masss flow of filtered dust";
    input Modelica.SIunits.Power powerConsumption "Auxiliary power";
 end Outline;

model Summary
     extends ClaRa.Basics.Icons.RecordIcon;
     Outline outline;
     ClaRa.Basics.Records.FlangeGas  inlet;
     ClaRa.Basics.Records.FlangeGas  outlet;
end Summary;

//## P A R A M E T E R S #######################################################################################
//_____________defintion of medium used in cell__________________________________________________________
  inner parameter TILMedia.GasTypes.BaseGas               medium = simCenter.flueGasModel
    "Medium to be used in tubes"                                                                                       annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));
  parameter Real epsilon_r = 10 "Dielectric number of flueGas";

  parameter Real specific_powerConsumption(unit="W.h/m3") = 0.15
    "Specific power consumption"                                                              annotation (Dialog(group="Fundamental Definitions"));

  inner parameter ClaRa.Basics.Units.MassFlowRate m_flow_nom= 10
    "Nominal mass flow rates at inlet"                                                              annotation(Dialog(tab="General", group="Nominal Values"));
  inner parameter ClaRa.Basics.Units.Pressure p_nom=1e5 "Nominal pressure"                    annotation(Dialog(group="Nominal Values"));
  inner parameter ClaRa.Basics.Units.EnthalpyMassSpecific h_nom=1e5
    "Nominal specific enthalpy"                                                                  annotation(Dialog(group="Nominal Values"));

  inner parameter ClaRa.Basics.Choices.Init initType=ClaRa.Basics.Choices.Init.noInit
    "Type of initialisation"                                                                                   annotation(Dialog(tab="Initialisation", choicesAllMatching));
  parameter ClaRa.Basics.Units.Temperature T_start= 380
    "Start value of system Temperature"                                                     annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.Pressure p_start= 1.013e5
    "Start value of system pressure"                                                      annotation(Dialog(tab="Initialisation"));
  parameter ClaRa.Basics.Units.MassFraction xi_start[medium.nc-1]={0.01,0,0.25,0,0.7,0,0,0.04,0}
    "Start value of system mass fraction"                                                   annotation(Dialog(tab="Initialisation"));
  inner parameter Boolean useHomotopy=simCenter.useHomotopy
    "True, if homotopy method is used during initialisation"                                                         annotation(Dialog(tab="Initialisation"));

// Quantaties required for the calculation of the separationRate
  parameter ClaRa.Basics.Units.Area A_filter = 100 "Collector area of E-Filter"
                                                                                annotation(Dialog(group="Geometry"));
  parameter ClaRa.Basics.Units.Length d_plate = 0.2
    "Distance  Plate-to-Plate or Plate-to-Wire, repectivaly"                                annotation(Dialog(group="Geometry"));
  parameter ClaRa.Basics.Units.Length diameter_particle = 50e-6
    "Average diameter of ash particles"                                                            annotation(Dialog(group="Geometry"));
  final parameter Real A1 = 1.257 "Auxiliary Area"
                                                  annotation(Dialog(group="Geometry"));
  final parameter Real A2 = 0.4 "Auxiliary Area" annotation(Dialog(group="Geometry"));
  final parameter Real A3 = 0.55 "Auxiliary Area"
                                                 annotation(Dialog(group="Geometry"));

  parameter Boolean use_dynamicMassbalance = true
    "True if a dynamic mass balance shall be applied"                                               annotation(Evaluate=true, Dialog(tab="Expert Settings"));

  parameter Boolean contributeToCycleSummary = simCenter.contributeToCycleSummary
    "True if component shall contribute to automatic efficiency calculation"                annotation(Dialog(tab="Summary and Visualisation"));
  parameter Boolean showData=true
    "True if a data port containing p,T,h,s,m_flow shall be shown, else false"              annotation (Dialog(tab="Summary and Visualisation"));
//## V A R I A B L E   P A R T##################################################################################

Real separationRate;
Modelica.SIunits.Velocity w_m
    "Migration speed of dust particles in the E-field";                          // Quantities required for the calculation of the particles' migration speed
Modelica.SIunits.Length lambda "Mean free path of particles";
Modelica.SIunits.ElectricFieldStrength E_applied
    "E-Field in Filter estimated as E = U/d with U being the applied potential and d the distance between elektrodes, refer to Riehle 1997";
     Modelica.Blocks.Interfaces.RealInput U_applied
    "Applied Voltage lading to E_applied=U_applied/d" annotation (Placement(
        transformation(
        extent={{-12,-12},{12,12}},
        rotation=270,
        origin={0,112}), iconTransformation(
        extent={{-12,-12},{12,12}},
        rotation=270,
        origin={-74,112})));

ClaRa.Basics.Units.DynamicViscosity mu_flueGas
    "Dynamic viscosity of flueGas in E-Filter";
Real Cu
    "Cunningham slip correction factor Cu = 1 + 2lambda/d *(A1 +A2exp[-A3*diameter_particle/lambda]])";

Modelica.SIunits.ElectricCharge Q_sat "Saturation charge of particles";
// Quantaties required for the calculation of the particles' saturation charge

ClaRa.Basics.Units.VolumeFlowRate V_flow
    "Volumeflow rate of flue Gas entering the E-Filter";

protected
   ClaRa.Basics.Units.EnthalpyMassSpecific h_out "Specific enthalpy at outlet";
   ClaRa.Basics.Units.EnthalpyMassSpecific h_in "Specific enthalpy at inlet";
   ClaRa.Basics.Units.EnthalpyMassSpecific h_dust
    "Specific enthalpy of separated dust";
   inner ClaRa.Basics.Units.EnthalpyMassSpecific h(start=TILMedia.GasFunctions.specificEnthalpy_pTxi(simCenter.flueGasModel, p_start, T_start, xi_start))
    "Specific enthalpy of gas";
   Real drhodt "Density derivative";
   Modelica.SIunits.Mass mass "Mass in component";
   Modelica.SIunits.Pressure p(start=p_start) "Pressure in component";
   Modelica.SIunits.MassFraction xi[medium.nc-1]( start=xi_start)
    "Mass fraction";
   Modelica.SIunits.MassFlowRate m_flow_dust_out "Mass flow of separated dust";
   Modelica.SIunits.Power powerConsumption "Power consumption";

//____Connectors________________________________________________________________________________________________
public
  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium) "Inlet port"
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortOut outlet(Medium=medium) "Outlet port"
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));
  ClaRa.Basics.Interfaces.HeatPort_a
                                   heat
    annotation (Placement(transformation(extent={{-10,90},{10,110}})));

//____replaceable models for heat transfer, pressure loss and geometry____________________________________________
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

  replaceable model Geometry =
      ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry
    constrainedby
    ClaRa.Basics.ControlVolumes.Fundamentals.Geometry.GenericGeometry
    "1st: choose geometry definition | 2nd: edit corresponding record"
    annotation (Dialog(group="Geometry"), choicesAllMatching=true);

public
  HeatTransfer heattransfer(heatSurfaceAlloc=1)
    annotation (Placement(transformation(extent={{-80,60},{-60,80}})));
  inner Geometry geo annotation (Placement(transformation(extent={{-48,60},
            {-28,80}})));

  PressureLoss pressureLoss annotation (Placement(transformation(extent={{12,60},
            {32,80}})));

//_____________________Media Objects_________________________________
protected
  TILMedia.Gas_pT     flueGasInlet(p=inlet.p,
  T=noEvent(actualStream(inlet.T_outflow)),
  xi=noEvent(actualStream(inlet.xi_outflow)),
  gasType = medium)
    annotation (Placement(transformation(extent={{-80,-20},{-60,0}})));
  TILMedia.Gas_pT     flueGasOutlet(p=outlet.p,
  T=noEvent(actualStream(outlet.T_outflow)),
  xi=noEvent(actualStream(outlet.xi_outflow)),
  gasType = medium)
    annotation (Placement(transformation(extent={{60,-20},{80,0}})));

  inner TILMedia.Gas_ph     bulk(
    computeTransportProperties=false,
    gasType = medium,p=p,h=h,xi=xi,
    stateSelectPreferForInputs=true)
    annotation (Placement(transformation(extent={{-10,-20},{10,0}})));
public
   TILMedia.Gas_pT     dust(
    computeTransportProperties=false,
    gasType = medium,
    stateSelectPreferForInputs=true,
    p=bulk.p,
    T=bulk.T,
    xi(start={if i==1 then 0.99999 else if i==5 then 0.00001 else 0 for i in 1:medium.nc-1})={if i==1  then 0.99999 else if i==5 then 0.00001 else 0 for i in 1:medium.nc-1})
    annotation (Placement(transformation(extent={{-10,-60},{10,-40}})));

Summary summary(outline(
    V=geo.volume,
    m=mass,
    H=mass*h,
    h=h,
    p=bulk.p,
    Delta_p=pressureLoss.Delta_p,
    lambda=lambda,
    Cu=Cu,
    w_m=w_m,
    E_applied=E_applied,
    Q_sat=Q_sat,
    separationRate=separationRate,
    m_flow_dust_out=m_flow_dust_out,
    powerConsumption=powerConsumption),
    inlet(m_flow=inlet.m_flow,
    T=inStream(inlet.T_outflow),
    p=inlet.p,
    h=flueGasInlet.h,
    xi = inStream(inlet.xi_outflow),
    H_flow = inlet.m_flow*flueGasInlet.h),
    outlet(m_flow=-outlet.m_flow,
    T=outlet.T_outflow,
    p=outlet.p,
    h=flueGasOutlet.h,
    xi = outlet.xi_outflow,
    H_flow = -outlet.m_flow*flueGasOutlet.h));

protected
  inner ClaRa.Basics.Records.IComGas_L2 iCom(
    m_flow_nom=m_flow_nom,
    T_bulk=bulk.T,
    p_bulk=bulk.p,
    fluidPointer_in=flueGasInlet.gasPointer,
    fluidPointer_bulk=bulk.gasPointer,
    fluidPointer_out=flueGasOutlet.gasPointer,
    mediumModel=medium,
    p_in=inlet.p,
    T_in=flueGasInlet.T,
    m_flow_in=inlet.m_flow,
    V_flow_in=0,
    xi_in=xi,
    p_out=outlet.p,
    T_out=flueGasOutlet.T,
    m_flow_out=outlet.m_flow,
    V_flow_out=0,
    xi_out=xi) annotation (Placement(transformation(extent={{80,-102},{100,-82}})));

public
  ClaRa.Basics.Interfaces.EyeOut eyeOut if showData annotation (
      Placement(transformation(extent={{80,-78},{120,-42}}),
        iconTransformation(extent={{90,-50},{110,-30}})));
protected
  ClaRa.Basics.Interfaces.EyeIn eye_int annotation (Placement(
        transformation(extent={{48,-68},{32,-52}}),
        iconTransformation(extent={{90,-84},{84,-78}})));
equation
  // Asserts ~~~~~~~~~~~~~~~~~~~
   assert(geo.volume>0, "The system volume must be greater than zero!");
    assert(geo.A_heat[heattransfer.heatSurfaceAlloc]>=0, "The area of heat transfer must be greater than zero!");

  // Port connection
  inlet.T_outflow  = bulk.T;
  outlet.T_outflow = bulk.T;

  inlet.xi_outflow  = xi;
  outlet.xi_outflow = xi;

  h_in=flueGasInlet.h;
  h_out=flueGasOutlet.h;
  //-----------------------------------------------------------
  h_dust = dust.h;
  mass = geo.volume * bulk.d;

  /*
  inlet.p=if useHomotopy then homotopy(p+pressureLoss.Delta_p + (geo.z_in-geo.z_out)*Modelica.Constants.g_n*bulk.d,
                                           p+pressureLoss.Delta_p + (geo.z_in-geo.z_out)*Modelica.Constants.g_n*d_nom)
              else p+pressureLoss.Delta_p + (geo.z_in-geo.z_out)*Modelica.Constants.g_n*bulk.d;*/

     inlet.p =  p+pressureLoss.Delta_p;// + (geo.z_in-geo.z_out)*Modelica.Constants.g_n*bulk.d;
     outlet.p = p;

  // Mass balance
  if use_dynamicMassbalance then
    inlet.m_flow + outlet.m_flow + m_flow_dust_out =  drhodt*geo.volume;
    der(xi) =
     1/mass * (inlet.m_flow*(flueGasInlet.xi - xi) + outlet.m_flow*(flueGasOutlet.xi - xi) + m_flow_dust_out*(dust.xi-xi));
  else
    inlet.m_flow + outlet.m_flow + m_flow_dust_out =  0;
    zeros(medium.nc-1) =
      (inlet.m_flow*(flueGasInlet.xi - xi) + outlet.m_flow*(flueGasOutlet.xi - xi) + m_flow_dust_out*(dust.xi-xi));
  end if;

  if inlet.m_flow > 0 and outlet.m_flow <=0 then
    V_flow = inlet.m_flow/flueGasInlet.d;
    m_flow_dust_out = separationRate *(-flueGasInlet.xi[1]*inlet.m_flow);
  else if  inlet.m_flow > 0 and outlet.m_flow > 0 then
        V_flow = (inlet.m_flow/flueGasInlet.d + outlet.m_flow/flueGasOutlet.d);
        m_flow_dust_out = separationRate * (-flueGasOutlet.xi[1]*outlet.m_flow-flueGasInlet.xi[1]*inlet.m_flow);
  else if inlet.m_flow <= 0 and outlet.m_flow <= 0 then
       V_flow = 1e-20;
       m_flow_dust_out = -1e-20;
       else
       V_flow = outlet.m_flow/flueGasOutlet.d;
       m_flow_dust_out = separationRate *(-flueGasOutlet.xi[1]*outlet.m_flow);
            end if;
      end if;
  end if;

 E_applied = U_applied/d_plate;

 // calculation of the migration velocity refering to Riehle "Basic and theoretical operation of ESPs" (1997)
 w_m = Q_sat*E_applied*Cu/(3.0*Modelica.Constants.pi*mu_flueGas*diameter_particle);
// with
 Q_sat = ((1.0+2.0*lambda/diameter_particle)^2 + (2.0/(1.0+2.0*lambda/diameter_particle))*(epsilon_r-1.0)/(epsilon_r+2.0)) * Modelica.Constants.pi*Modelica.Constants.epsilon_0*diameter_particle^2*E_applied;
//and
 Cu = 1 + 2.0*lambda/diameter_particle*(A1 + A2*Modelica.Math.exp(-A3*diameter_particle/lambda));
 //with lambda the mean free path of molecules/paticles in the flueGas (check, which diameter is physical plausible)
 lambda = 1.0/(4.0*Modelica.Constants.pi*sqrt(2.0)*diameter_particle^2*Modelica.Constants.N_A/(Modelica.Constants.R*flueGasInlet.T/flueGasInlet.p));

 mu_flueGas  = flueGasInlet.transp.eta;

  separationRate =  1- Modelica.Math.exp(-w_m*A_filter/V_flow);
  powerConsumption = V_flow * specific_powerConsumption*3600.;

  if use_dynamicMassbalance then
    drhodt = bulk.drhodh_pxi * der(h)
             + bulk.drhodp_hxi * der(p)
             + sum({bulk.drhodxi_ph[i] * der(bulk.xi[i]) for i in 1:medium.nc-1});
  else
     drhodt = bulk.drhodh_pxi * der(h)
             + bulk.drhodp_hxi * der(p);
  end if;

       der(h) =  (inlet.m_flow*(h_in-h) + outlet.m_flow*(h_out-h) + m_flow_dust_out*(h_dust-h)  + geo.volume*der(p) + heat.Q_flow)/mass;

  //______________Eye port variable definition________________________
  eye_int.m_flow = -outlet.m_flow;
  eye_int.T = bulk.T-273.15;
  eye_int.s = bulk.s/1e3;
  eye_int.p = bulk.p/1e5;
  eye_int.h = bulk.h/1e3;

initial equation

  if initType == ClaRa.Basics.Choices.Init.steadyState then
      der(h)=0;
      der(p)=0;
  elseif initType == ClaRa.Basics.Choices.Init.steadyPressure then
      der(p)=0;
  elseif initType == ClaRa.Basics.Choices.Init.steadyEnthalpy then
      der(h)=0;
  elseif initType == ClaRa.Basics.Choices.Init.noInit then
      bulk.T=T_start;
    end if;

equation
  connect(heattransfer.heat, heat) annotation (Line(
      points={{-60,70},{-54,70},{-54,90},{0,90},{0,100}},
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
            extent={{-100,-100},{100,100}}),
                                        graphics),
    Documentation(info="<html>
<p><b>Model description: </b>An detailed E-filter model</p>
<p><b>Contact:</b> Andre Th&uuml;ring, Lasse Nielsen, TLK-Thermo GmbH</p>
<p><b>FEATURES</b> </p>
<p><ul>
<li>This model uses TILMedia</li>
<li>Calculates the separation rate with migration speed and a mean free path for molecules and the Cunningham&nbsp;slip&nbsp;correction&nbsp;factor&nbsp;Cu&nbsp;=&nbsp;1&nbsp;+&nbsp;2lambda/d&nbsp;*(A1&nbsp;+A2exp[-A3*diameter_particle/lambda]])</li>
<li>Calculation&nbsp;of&nbsp;the&nbsp;migration&nbsp;velocity&nbsp;refering&nbsp;to&nbsp;Riehle&nbsp;&QUOT;Basic&nbsp;and&nbsp;theoretical&nbsp;operation&nbsp;of&nbsp;ESPs&QUOT;&nbsp;(1997)</li>
<li>Power consumption is calculated with a given specific power consumption </li>
<li>Stationary mass and energy balance</li>
</ul></p>
</html>"));
end E_Filter_L2_detailed;
