within Exergy.XClaRa.Components.TurboMachines.Compressors.Fundamentals;
model GetInputsHydraulic "Get enabled inputs and parameters of disabled inputs"
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

  extends Modelica.Blocks.Interfaces.BlockIcon;

  Modelica.Blocks.Interfaces.RealInput m_flow_in "Prescribed mass flow rate"
    annotation (Placement(transformation(extent={{-140,-100},{-100,-60}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealInput V_flow_in "Prescribed volume flow rate"
    annotation (Placement(transformation(extent={{-140,-40},{-100,0}},
          rotation=0)));
  Modelica.Blocks.Interfaces.RealInput dp_in "Prescribed pressure increase"
    annotation (Placement(transformation(extent={{-140,60},{-100,100}},rotation=
           0)));
  Modelica.Blocks.Interfaces.RealInput P_shaft_in "Power input"
                                         annotation (Placement(transformation(
        extent={{20,20},{-20,-20}},
        rotation=180,
        origin={-120,20})));

  annotation (Diagram(graphics));
end GetInputsHydraulic;
