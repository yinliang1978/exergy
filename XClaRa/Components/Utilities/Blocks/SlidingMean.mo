within Exergy.XClaRa.Components.Utilities.Blocks;
block SlidingMean "Sliding  mean value"
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                            //
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

  extends Modelica.Blocks.Interfaces.SISO;

  parameter ClaRa.Basics.Units.Time Tau_mean= 1 "Time interval for averaging";

protected
  ClaRa.Basics.Functions.ClaRaDelay.ExternalTable pointer_Int = ClaRa.Basics.Functions.ClaRaDelay.ExternalTable();
  Real integral "Integral of u over time at time";
  Real integral_old(start=0) "Integral of at time - Tau_mean";
  parameter Real u_start(fixed=false) "Assumed u between time=[-Tau_mean, 0]";

equation
  der(integral) = u;
  integral_old = ClaRa.Basics.Functions.ClaRaDelay.getDelayValuesAtTime(pointer_Int, time, integral,  time-Tau_mean) - (min(1,max(0,(Tau_mean-time)/Tau_mean)))*u_start*Tau_mean;
  y=(integral-integral_old)/Tau_mean;

initial equation
  integral=0;
  u_start=u;
  annotation (Icon(graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(
          points={{-78,0},{-56,-4},{-30,0},{-4,-10},{26,4},{76,-4}},
          color={27,36,42},
          smooth=Smooth.Bezier),
        Line(
          points={{-78,0},{-60,-20},{-42,30},{-32,-14},{-20,-14},{-4,-64},{4,24},{8,-28},{24,76},{36,-16},{38,-54},{69.9688,33.9141},{72,-12}},
          color={221,222,223},
          smooth=Smooth.Bezier),
        Line(
          points={{-78,80},{-78,-78},{78,-78}},
          color={221,222,223})}));
end SlidingMean;
