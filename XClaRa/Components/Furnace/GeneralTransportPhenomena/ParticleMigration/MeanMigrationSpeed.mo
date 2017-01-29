within Exergy.XClaRa.Components.Furnace.GeneralTransportPhenomena.ParticleMigration;
model MeanMigrationSpeed
  "Determines the mean migration speed in dependence of the fluegas flow rates"
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
    ClaRa.Components.Furnace.GeneralTransportPhenomena.ParticleMigration.PartialMigrationSpeed;

    parameter ClaRa.Basics.Units.Velocity
                                      w_initial = 2 annotation(Dialog(tab="Start and initialisation values"));
    outer ClaRa.Basics.Units.Time Tau;
initial equation

  w = w_initial;

equation
 der(w) = 1/Tau*((V_flow_flueGas_in-V_flow_flueGas_out)/2.0/A_cross - w);

  annotation (Documentation(info="<html>
<p><b>Model description: </b>Model for a mean migration speed calculation according to flue gas volume flow rates</p>
<p><b>Contact:</b> Lasse Nielsen, TLK-Thermo GmbH</p>
</html>"));
end MeanMigrationSpeed;
