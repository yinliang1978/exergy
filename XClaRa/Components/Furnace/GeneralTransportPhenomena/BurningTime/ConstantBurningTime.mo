within Exergy.XClaRa.Components.Furnace.GeneralTransportPhenomena.BurningTime;
model ConstantBurningTime "Assumes a constant burning time"
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

extends
    ClaRa.Components.Furnace.GeneralTransportPhenomena.BurningTime.PartialBurningTime;
    parameter ClaRa.Basics.Units.Time
                                   Tau_burn_const = 1.5;

equation
t = Tau_burn_const;
  annotation (Documentation(info="<html>
<p><b>Model description: </b>Model for constant burning time</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));
end ConstantBurningTime;
