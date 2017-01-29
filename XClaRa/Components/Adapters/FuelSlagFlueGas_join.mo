within Exergy.XClaRa.Components.Adapters;
model FuelSlagFlueGas_join
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
  extends ClaRa.Basics.Icons.Adapter3_bw;
 //__________________________/ Media definintions \______________________________________________
  outer ClaRa.SimCenter simCenter;
  inner parameter ClaRa.Basics.Media.Fuel.PartialFuel fuelType=simCenter.fuelModel1
    "Fuel elemental composition used for combustion" annotation(choicesAllMatching, Dialog(group="Fundamental Medium Definitions"));
  inner parameter ClaRa.Basics.Media.Fuel.PartialSlag slagType=simCenter.slagModel
    "Slag properties" annotation(choicesAllMatching, Dialog(group="Fundamental Medium Definitions"));
  inner parameter TILMedia.GasTypes.BaseGas               flueGas = simCenter.flueGasModel
    "Medium to be used in tubes"  annotation(choicesAllMatching, Dialog(group="Fundamental Medium Definitions"));

  ClaRa.Basics.Interfaces.Fuel_inlet fuel_inlet(fuelType=fuelType)
    annotation (Placement(transformation(extent={{-110,50},{-90,70}})));
  ClaRa.Basics.Interfaces.GasPortIn flueGas_inlet(Medium=flueGas)
    annotation (Placement(transformation(extent={{-110,-70},{-90,-50}})));
  ClaRa.Basics.Interfaces.Slag_outlet slag_outlet(slagType=slagType)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}}),
        iconTransformation(extent={{-110,-10},{-90,10}})));

  ClaRa.Basics.Interfaces.FuelSlagFlueGas_outlet fuelSlagFlueGas_outlet(
    flueGas(Medium=flueGas),
    final fuelType=fuelType,
    final slagType=slagType) annotation (Placement(transformation(
        extent={{-10,-10},{10,10}},
        rotation=270,
        origin={100,0}), iconTransformation(
        extent={{-10,-10},{10,10}},
        rotation=0,
        origin={100,0})));

equation
  fuelSlagFlueGas_outlet.flueGas.m_flow = -flueGas_inlet.m_flow;
  fuelSlagFlueGas_outlet.flueGas.T_outflow = inStream(flueGas_inlet.T_outflow);
  flueGas_inlet.T_outflow = inStream(fuelSlagFlueGas_outlet.flueGas.T_outflow);
  fuelSlagFlueGas_outlet.flueGas.xi_outflow = inStream(flueGas_inlet.xi_outflow);
  flueGas_inlet.xi_outflow = inStream(fuelSlagFlueGas_outlet.flueGas.xi_outflow);
  fuelSlagFlueGas_outlet.flueGas.p = flueGas_inlet.p;

  fuelSlagFlueGas_outlet.fuel.m_flow = -fuel_inlet.m_flow;
  fuelSlagFlueGas_outlet.fuel.T_outflow = inStream(fuel_inlet.T_outflow);
  fuel_inlet.T_outflow = inStream(fuelSlagFlueGas_outlet.fuel.T_outflow);
  fuelSlagFlueGas_outlet.fuel.xi_outflow = inStream(fuel_inlet.xi_outflow);
  fuel_inlet.xi_outflow = inStream(fuelSlagFlueGas_outlet.fuel.xi_outflow);
  fuelSlagFlueGas_outlet.fuel.LHV_outflow = inStream(fuel_inlet.LHV_outflow);
  fuel_inlet.LHV_outflow = inStream(fuelSlagFlueGas_outlet.fuel.LHV_outflow);
  fuelSlagFlueGas_outlet.fuel.cp_outflow = inStream(fuel_inlet.cp_outflow);
  fuel_inlet.cp_outflow = inStream(fuelSlagFlueGas_outlet.fuel.cp_outflow);
  fuelSlagFlueGas_outlet.fuel.p = fuel_inlet.p;
  fuelSlagFlueGas_outlet.fuel.LHV_calculationType= fuel_inlet.LHV_calculationType;

  fuelSlagFlueGas_outlet.slag.m_flow = -slag_outlet.m_flow;
  fuelSlagFlueGas_outlet.slag.T_outflow = inStream(slag_outlet.T_outflow);
  slag_outlet.T_outflow = inStream(fuelSlagFlueGas_outlet.slag.T_outflow);
  fuelSlagFlueGas_outlet.slag.p = slag_outlet.p;

  annotation (Icon(graphics),               Diagram(graphics));
end FuelSlagFlueGas_join;
