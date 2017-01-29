within Exergy.XClaRa.Components.BoundaryConditions.Check;
model TestFlueGasComposition
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
  extends ClaRa.Basics.Icons.PackageIcons.ExecutableExampleb60;
  GasCompositionByMassFractions flueGasComposition(
    xi_ASH=0.05,
    xi_CO=0.05,
    xi_CO2=0.5,
    xi_SO2=0.0,
    xi_N2=0.1,
    xi_O2=0.1,
    xi_NO=0.0,
    xi_H2O=0.2,
    xi_NH3=0)
    annotation (Placement(transformation(extent={{-100,14},{-80,34}})));
  BoundaryGas_Txim_flow idealGasFlowSource_XRG(
    variable_xi=true,
    medium=simCenter.flueGasModel,
    T_const=650 + 273.15,
    m_flow_const=1) annotation (Placement(transformation(extent={{-60,20},{-40,40}})));
  BoundaryGas_pTxi idealGasPressureSink_XRG(
    medium=simCenter.flueGasModel,
    variable_p=false,
    variable_xi=false,
    p_const=100000) annotation (Placement(transformation(extent={{80,20},{60,40}})));
  inner ClaRa.SimCenter simCenter annotation (Placement(
        transformation(extent={{80,80},{100,100}})));
equation
  connect(flueGasComposition.X, idealGasFlowSource_XRG.xi) annotation (Line(
      points={{-78,24},{-60,24}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(idealGasFlowSource_XRG.gas_a, idealGasPressureSink_XRG.gas_a)
    annotation (Line(
      points={{-40,30},{60,30}},
      color={118,106,98},
      thickness=0.5,
      smooth=Smooth.None));
  annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                      graphics));
end TestFlueGasComposition;
