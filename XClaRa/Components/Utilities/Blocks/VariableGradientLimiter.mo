within Exergy.XClaRa.Components.Utilities.Blocks;
block VariableGradientLimiter
  "Limit the range of a signal with variable limits"
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
  import ClaRa.Basics.Functions.Stepsmoother;
  parameter Real Nd(min=Modelica.Constants.small) = 10
    "|Expert Settings|Input - Output Coupling|The higher Nd, the closer y follows u";

  parameter Boolean useThresh=false
    "|Expert Settings|Numerical Noise Suppression|Use threshould for suppression of numerical noise";
  parameter Real thres(max=1e-6)=1e-7 "If abs(u-y)< thres, y becomes a simple pass through of u. Increasing thres can improve simulation speed. However to large values can make the simulation unstable. 
     A good starting point is the choice thres = tolerance/1000."                         annotation (Dialog(enable = useThresh,tab="Expert Settings",group="Numerical Noise Suppression"));
  parameter Boolean constantLimits= false
    "True, if gradient limits are constant";
  parameter Real maxGrad_const = 1 "Constant max gradient" annotation (Dialog(enable = constantLimits));
  parameter Real minGrad_const = -1 "Constant min gradient" annotation (Dialog(enable = constantLimits));

  Modelica.Blocks.Interfaces.RealInput maxGrad(value=maxGrad_) if not constantLimits
    "Maximum Gradient allowd" annotation (Placement(transformation(extent={{
            -140,60},{-100,100}}, rotation=0)));
  Modelica.Blocks.Interfaces.RealInput minGrad(value=minGrad_) if not constantLimits
    "Minimum Gradient allowd" annotation (Placement(transformation(extent={{
            -140,-100},{-100,-60}}, rotation=0)));

protected
  Real maxGrad_;
  Real minGrad_;
  Real y_aux;
equation

  if constantLimits then
    maxGrad_ = maxGrad_const;
    minGrad_ = minGrad_const;
  end if;
  der(y_aux) = smooth(1,noEvent(min(maxGrad_,max(minGrad_,(u-y_aux)*Nd))));

  if useThresh then
     y= Stepsmoother(1, 0.1, abs(der(y_aux))/thres)*y_aux + (1-Stepsmoother(1, 0.1, abs(der(y_aux))/thres))*u;
  else
    y=y_aux;
  end if;

initial equation
  y_aux=u;

  annotation (
    Documentation(info="<html>
<p>
The Limiter block passes its input signal as output signal
as long as the input is within the upper and lower
limits specified by the two additional inputs limit1 and
limit2. If this is not the case, the corresponding limit
is passed as output.
</p>
</HTML>
"), Icon(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{100,100}}),
                    graphics={
        Rectangle(
          extent={{-100,100},{100,-100}},
          lineColor={221,222,223},
          fillColor={118,124,127},
          fillPattern=FillPattern.Solid),
        Line(points={{-80,-80},{-80,82}},
                                      color={221,222,223}),
        Line(points={{-80,-80},{80,-80}},
                                      color={221,222,223}),
        Text(
          extent={{-150,150},{150,110}},
          lineColor={27,36,42},
          textString="%name"),
        Line(points={{-80,-60},{-18,-4},{80,-60}},
                                                 color={27,36,42}),
        Line(points={{-80,-60},{0,56},{80,-60}}, color={221,222,223})}),
    Diagram(coordinateSystem(preserveAspectRatio=true, extent={{-100,-100},{
            100,100}})));
end VariableGradientLimiter;
