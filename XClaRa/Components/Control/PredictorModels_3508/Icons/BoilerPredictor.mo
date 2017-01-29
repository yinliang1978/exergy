within Exergy.XClaRa.Components.Control.PredictorModels_3508.Icons;
model BoilerPredictor "An icon for a boiler model - relative units ('p.u.')"
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

  annotation (Icon(graphics={Polygon(
          points={{-4,-48},{28,-100},{72,-100},{100,-48},{100,100},{-16,100},{
              -62,66},{-40,38},{-4,68},{-4,-48}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={135,135,135},
          fillPattern=FillPattern.Solid), Polygon(
          points={{-100,-10},{-12,-10},{-38,10},{-32,10},{0,-14},{-32,-38},{-40,
              -38},{-10,-16},{-100,-16},{-100,-10}},
          lineColor={0,0,0},
          smooth=Smooth.None,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Polygon(
          points={{0,-16},{26,-30},{54,-26},{74,-12},{74,18},{54,36},{54,50},{
              50,58},{48,46},{40,32},{40,14},{32,4},{20,-4},{6,-10},{0,-16}},
          lineColor={0,0,0},
          smooth=Smooth.Bezier,
          fillColor={227,227,227},
          fillPattern=FillPattern.Solid)}));
end BoilerPredictor;
