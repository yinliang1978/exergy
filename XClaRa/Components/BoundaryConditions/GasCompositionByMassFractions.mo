within Exergy.XClaRa.Components.BoundaryConditions;
model GasCompositionByMassFractions "set (flue) gas composition graphically"
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

  extends ClaRa.Basics.Icons.MassComposition;
  Modelica.Blocks.Interfaces.RealOutput X[medium.nc - 1]
    "composition of gas to be set"
    annotation (Placement(transformation(extent={{120,0},{140,20}}),
        iconTransformation(extent={{100,-20},{140,20}})));
  parameter Real xi_ASH "Mass fraction of ash";
  parameter Real xi_CO "Mass fraction of carbon monoxide (CO)";
  parameter Real xi_CO2 "Mass fraction of carbon dioxide (CO2)";
  parameter Real xi_SO2 "Mass fraction of sulphur dioxide (SO2)";
  parameter Real xi_N2 "Mass fraction of nitrogen (N2)";
  parameter Real xi_O2 "Mass fraction of oxygen (O2)";
  parameter Real xi_NO "Mass fraction of nitrogen oxide (NO)";
  parameter Real xi_H2O "Mass fraction of water (H2O)";
  parameter Real xi_NH3 "Mass fraction of ammonia (NH3)";
  Real sumXi "control sum of set mass fractions";

protected
  outer ClaRa.SimCenter simCenter;
public
  Modelica.SIunits.MassFraction xi_in[medium.nc - 1];
  TILMedia.GasTypes.BaseGas      medium = simCenter.flueGasModel;
equation
  xi_in[1] = xi_ASH;
  xi_in[2] = xi_CO;
  xi_in[3] = xi_CO2;
  xi_in[4] = xi_SO2;
  xi_in[5] = xi_N2;
  xi_in[6] = xi_O2;
  xi_in[7] = xi_NO;
  xi_in[8] = xi_H2O;
  xi_in[9] = xi_NH3;

  X = xi_in;
algorithm
  for i in 1:  1:  size(xi_in,1) loop
    sumXi :=sumXi + xi_in[i];
  end for;

annotation (Documentation(info="<html>
<p>
<a href=\"modelica://ClaRa/figures/GasCompositionByMassFractions.html\" title=\"modelica://ClaRa/figures/GasCompositionByMassFractions.html\">
GasCompositionByMassFractions Model</a><br>
</html>"),                                                        Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},{100,100}}),
                                                                       graphics),
                                   Diagram(graphics),
              Icon(graphics));
end GasCompositionByMassFractions;
