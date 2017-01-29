within Exergy.XClaRa.Components.Utilities.Blocks.Fundamentals;
function random
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
  extends ClaRa.Basics.Icons.Function;

  input Real si[3] "Input random seed";
  output Real x "Uniform random variate between 0 and 1";
  output Real so[3] "Output random seed";

algorithm
  so[1] := abs(rem((171 * si[1]),30269));
  so[2] := abs(rem((172 * si[2]),30307));
  so[3] := abs(rem((170 * si[3]),30323));

// zero is a poor Seed, therfore substitute 1;

  if so[1] == 0 then
    so[1] := 1;
  end if;

  if so[2] == 0 then
    so[2] := 1;
  end if;

  if so[3] == 0 then
    so[3] := 1;
  end if;

  x := rem((so[1]/30269.0 +so[2]/30307.0 +
  so[3]/30323.0),1.0);

end random;
