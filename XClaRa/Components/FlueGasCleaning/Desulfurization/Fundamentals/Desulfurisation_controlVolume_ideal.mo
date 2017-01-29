within Exergy.XClaRa.Components.FlueGasCleaning.Desulfurization.Fundamentals;
model Desulfurisation_controlVolume_ideal
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
    extends ClaRa.Basics.Icons.Box;

  outer ClaRa.SimCenter simCenter;

  inner parameter TILMedia.GasTypes.BaseGas  medium = simCenter.flueGasModel
    "Medium to be used in tubes" annotation(choicesAllMatching, Dialog(group="Fundamental Definitions"));

  ClaRa.Basics.Interfaces.GasPortIn inlet(Medium=medium)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));
  ClaRa.Basics.Interfaces.GasPortOut outlet(Medium=medium)
    annotation (Placement(transformation(extent={{90,-10},{110,10}}),
        iconTransformation(extent={{90,-10},{110,10}})));

//  final parameter Modelica.SIunits.MolarInternalEnergy delta_f_H_CaCO3 = -1207.1e3;
//  final parameter Modelica.SIunits.MolarInternalEnergy delta_f_H_CaSO4_H2O = -2023e3;

final parameter Modelica.SIunits.MolarMass M_CaSO4_H2O=0.172141
    "Molar mass of gypsum";
final parameter Modelica.SIunits.MolarMass M_CaCO3=0.10009
    "Molar mass of calcium carbonate";

Real m_flow_aux;
parameter Boolean useStabilisedMassFlow=false
    "|Expert Settings|Numerical Robustness|";
    parameter SI.Time Tau= 0.001 "Time Constant of Stabilisation" annotation(Dialog(tab="Expert Settings", group = "Numerical Robustness", enable=useStabilisedMassFlow));

parameter Real SOx_separationRate = 0.95 "Efficiency of SOx separation";
parameter Modelica.SIunits.Temperature T_in_H2O = 313.15
    "Inlet Temperature of water";
parameter Modelica.SIunits.MassFraction xi_start[medium.nc-1]=zeros(medium.nc-1)
    "Start value of system mass fraction" annotation(Dialog(tab="Initialisation"));

//required molar flow rates of reaction educts
Modelica.SIunits.MolarFlowRate n_flow_CaCO3_req
    "Required molar flow of calcium carbonate";
Modelica.SIunits.MolarFlowRate n_flow_O2_req
    "Additional required molar flow of oxygen";
Modelica.SIunits.MolarFlowRate n_flow_H2O_req "Required molar flow of water";

//molar flow rates of reaction educts inside flue gas
Modelica.SIunits.MolarFlowRate n_flow_SO2_in
    "Molar flow rate of sulfur dioxide at inlet";
Modelica.SIunits.MolarFlowRate n_flow_O2_in
    "Molar flow rate of oxygen at inlet";
Modelica.SIunits.MolarFlowRate n_flow_H2O_in
    "Molar flow rate of water at inlet";

//molar flow rates of products
Modelica.SIunits.MolarFlowRate n_flow_CaSO4_H2O_out
    "Molar flow rate of gypsum outlet (no connector)";
Modelica.SIunits.MolarFlowRate n_flow_CO2_out
    "Molar flow rate of carbon dioxide at outlet";
Modelica.SIunits.MolarFlowRate n_flow_H2O_out(start=1)
    "Molar flow rate of water at outlet";
Modelica.SIunits.MolarFlowRate n_flow_H2O_sep
    "Molar flow rate of separated water (no connector)";

Modelica.SIunits.MassFlowRate m_flow_SOx_sep
    "Mass flow of separated sulfur dioxide";
Modelica.SIunits.MassFlowRate m_flow_CaSO4_H2O_out
    "Mass flow of gypsum (no connector)";
Modelica.SIunits.MassFlowRate m_flow_H2O_req "Mass flow of required water";
Modelica.SIunits.MassFlowRate m_flow_O2_req "Mass flow of required oxygen";
Modelica.SIunits.MassFlowRate m_flow_CaCO3_req
    "Mass flow of required calcium carbonate";
Modelica.SIunits.MassFlowRate m_flow_O2_sep "Mass flow of separated oxygen";
Modelica.SIunits.MassFlowRate m_flow_H2O_sep "Mass flow of separated water";
Modelica.SIunits.MassFlowRate m_flow_CO2_prod
    "Mass flow of produced carbon dioxide";

ClaRa.Basics.Units.EnthalpyMassSpecific h_out "Specific enthalpy at outlet";
//ClaRa.Basics.Units.EnthalpyMassSpecific h_out_del "Pseudo state for specific enthalpy at outlet";

     TILMedia.GasObjectFunctions.GasPointer GasPointerOutlet=
      TILMedia.GasObjectFunctions.GasPointer(
       medium.concatGasName,
       8,
       medium.xi_default,
       medium.nc_propertyCalculation,
       medium.nc,
       0,
       0);

  TILMedia.Gas_pT     flueGasInlet(p=inlet.p,
  T=inStream(inlet.T_outflow),
  xi=inStream(inlet.xi_outflow),
  gasType = medium)
    annotation (Placement(transformation(extent={{-78,-12},{-58,8}})));

  TILMedia.Gas_ph     flueGasOutlet(
  gasType = medium,
    p=inlet.p,
    h=h_out,
    xi(start= xi_start)=outlet.xi_outflow)
    annotation (Placement(transformation(extent={{62,-12},{82,8}})));

initial equation
  - m_flow_aux= inlet.m_flow - m_flow_SOx_sep + m_flow_O2_req - m_flow_O2_sep + m_flow_H2O_req + m_flow_CO2_prod + m_flow_CaCO3_req - m_flow_CaSO4_H2O_out;
 // h_out_del = h_out;

equation
n_flow_SO2_in =inlet.m_flow*inStream(inlet.xi_outflow[4])/flueGasInlet.M_i[4];
n_flow_O2_in =inlet.m_flow*inStream(inlet.xi_outflow[6])/flueGasInlet.M_i[6];
n_flow_H2O_in =inlet.m_flow*inStream(inlet.xi_outflow[8])/flueGasInlet.M_i[8];
n_flow_H2O_out = - outlet.m_flow * flueGasOutlet.xi_s/flueGasInlet.M_i[8];
n_flow_H2O_sep = 2 * SOx_separationRate * n_flow_SO2_in;

n_flow_CaCO3_req = SOx_separationRate * n_flow_SO2_in;

n_flow_O2_req =
if n_flow_O2_in > 0.5 * SOx_separationRate * n_flow_SO2_in then
 0 else
     0.5 * SOx_separationRate * n_flow_SO2_in - n_flow_O2_in;

n_flow_H2O_req =
if n_flow_H2O_in > 2 * SOx_separationRate * n_flow_SO2_in + n_flow_H2O_out then
  0 else
  n_flow_H2O_sep + n_flow_H2O_out - n_flow_H2O_in;

n_flow_CaSO4_H2O_out = SOx_separationRate * n_flow_SO2_in;
n_flow_CO2_out = SOx_separationRate * n_flow_SO2_in;

m_flow_SOx_sep =SOx_separationRate*inlet.m_flow*inStream(inlet.xi_outflow[4]);
m_flow_O2_sep = 0.5 * SOx_separationRate * n_flow_SO2_in * flueGasInlet.M_i[6];
m_flow_H2O_sep = n_flow_H2O_sep * flueGasInlet.M_i[8];

m_flow_CaSO4_H2O_out = n_flow_CaSO4_H2O_out * M_CaSO4_H2O; //Inherits the separated H2O
m_flow_H2O_req = n_flow_H2O_req * flueGasInlet.M_i[8];
m_flow_O2_req = n_flow_O2_req *flueGasInlet.M_i[6];
m_flow_CaCO3_req = n_flow_CaCO3_req * M_CaCO3;
m_flow_CO2_prod = 1*SOx_separationRate*n_flow_SO2_in*flueGasInlet.M_i[3];

if (useStabilisedMassFlow==false) then
  - outlet.m_flow = inlet.m_flow - m_flow_SOx_sep + m_flow_O2_req - m_flow_O2_sep + m_flow_H2O_req + m_flow_CO2_prod + m_flow_CaCO3_req - m_flow_CaSO4_H2O_out; //Without H2O_sep because it is inherited within m_flow_CaSO4_H2O_out
  der(m_flow_aux) = 0;
else
  der(m_flow_aux) = 1/Tau *(-(inlet.m_flow - m_flow_SOx_sep + m_flow_O2_req - m_flow_O2_sep + m_flow_H2O_req + m_flow_CO2_prod + m_flow_CaCO3_req - m_flow_CaSO4_H2O_out) -m_flow_aux);
  outlet.m_flow = m_flow_aux;
end if;

//Komponentenbilanz
 for i in 1:(medium.nc-1) loop
    if i == 3 then
      inlet.m_flow*inStream(inlet.xi_outflow[3]) + (1*SOx_separationRate*
        n_flow_SO2_in*flueGasInlet.M_i[3]) = -outlet.m_flow*outlet.xi_outflow[3];
    else if i == 4 then
        inlet.m_flow*inStream(inlet.xi_outflow[4]) - (1*SOx_separationRate*
          n_flow_SO2_in*flueGasInlet.M_i[4]) = -outlet.m_flow*outlet.xi_outflow[
          4];
    else if i == 6 then
     if n_flow_O2_in < 0.5 * SOx_separationRate * n_flow_SO2_in then
            outlet.xi_outflow[6] = 0;
                                 else
            inlet.m_flow*inStream(inlet.xi_outflow[6]) - (0.5*
              SOx_separationRate*n_flow_SO2_in*flueGasInlet.M_i[6]) = -outlet.m_flow
              *outlet.xi_outflow[6];
     end if;
   else if i == 8 then
            outlet.xi_outflow[8] = flueGasOutlet.xi_s;
                                                  //Outlet flue gas is fully saturated with water
       //outlet.Xi_outflow[8] = xi_s_del;
   else     inlet.m_flow*inStream(inlet.xi_outflow[i]) = -outlet.m_flow*outlet.xi_outflow[
              i];
   end if;
   end if;
   end if;
   end if;
 end for;

0 = inlet.m_flow * flueGasInlet.h + (m_flow_H2O_req-m_flow_H2O_sep) * (-flueGasInlet.delta_hv) + outlet.m_flow * h_out;
//der(h_out_del)=1/Tau*(h_out-h_out_del);
outlet.T_outflow= flueGasOutlet.T;

inlet.p = outlet.p;

  inlet.xi_outflow
                 = inStream(outlet.xi_outflow);
inlet.T_outflow = inStream(outlet.T_outflow);
  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}),
                   graphics));
end Desulfurisation_controlVolume_ideal;
