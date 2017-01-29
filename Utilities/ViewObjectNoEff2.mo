within Exergy.Utilities;
model ViewObjectNoEff2

  import Exergy.Utilities.Types.*;

  //outer Exergy.Utilities.RefEnv refEnv "reference enviroment properties";
  //parameter Integer nHeat annotation (HideResult=true);
  parameter Integer[Exergy.Utilities.Types.EnergyTypes] nEnergy;

  Exergy.Utilities.EnergyRateData h[nEnergy[EnergyTypes.Enthalpy]];
  Exergy.Utilities.EnergyRateData e[nEnergy[EnergyTypes.InternalEnergy]];
  Exergy.Utilities.EnergyRateData w[nEnergy[EnergyTypes.Power]];
  Exergy.Utilities.EnergyRateData q[nEnergy[EnergyTypes.HeatTransfer]];

  Exergy.Utilities.EnergyData H[nEnergy[EnergyTypes.Enthalpy]];
  Exergy.Utilities.EnergyData E[nEnergy[EnergyTypes.InternalEnergy]];
  Exergy.Utilities.EnergyData W[nEnergy[EnergyTypes.Power]];
  Exergy.Utilities.EnergyData Q[nEnergy[EnergyTypes.HeatTransfer]];

  Modelica.SIunits.EnergyFlowRate E_flow_bal "Energy balance error";
  Modelica.SIunits.EnergyFlowRate Ex_flow_loss "exergy loss";
  Modelica.SIunits.InternalEnergy  Ex_loss "exergy loss";

  Exergy.Utilities.ViewPort viewOutput
    annotation (Placement(transformation(extent={{90,-10},{110,10}})));

//protected

equation
  h.E_flow=der(H.E);
  h.Ex_flow=der(H.Ex);
  e.E_flow=der(E.E);
  e.Ex_flow=der(E.Ex);
  w.E_flow=der(W.E);
  w.Ex_flow=der(W.Ex);
  q.E_flow=der(Q.E);
  q.Ex_flow=der(Q.Ex);

  Ex_flow_loss=der(Ex_loss);

  -E_flow_bal = sum(e.E_flow) - (sum(h.E_flow) + sum(w.E_flow) + sum(q.E_flow));
  -Ex_flow_loss = sum(e.Ex_flow) - (sum(h.Ex_flow) + sum(w.Ex_flow) + sum(q.Ex_flow));

 //-Ex_loss=sum(E.Ex) - (sum(H.Ex) + sum(W.Ex) + sum(Q.Ex));

  //volume properties summation
  viewOutput.E.E = sum(E.E);
  viewOutput.E.Ex = sum(E.Ex);

 // if noEvent(Ex_flow_input>0) then
    //  Ex_efficient*Ex_flow_input=Ex_loss;

 // else
   // Ex_efficient=0;
 // end if;

 // viewOutput.internalEnergy.E = der(viewOutput.internalEnergy.E_int);
 // viewOutput.internalEnergy.Ex = der(viewOutput.internalEnergy.Ex_int);

  annotation (
    defaultComponentName="viewObject",
    Icon(coordinateSystem(
        preserveAspectRatio=true,
        extent={{-100.0,-100.0},{100.0,100.0}},
        initialScale=0.1), graphics={Text(
          extent={{4,10},{-6,-16}},
          lineColor={28,108,200},
          fontSize=28,
          textString="*")}),
    Diagram(coordinateSystem(
        preserveAspectRatio=false,
        extent={{-100,-100},{100,100}},
        initialScale=0.1), graphics={Text(
          lineColor={0,0,127},
          extent={{30.0,60.0},{30.0,110.0}},
          textString="%name"), Text(
          extent={{0,2},{0,-4}},
          lineColor={28,108,200},
          fontSize=28,
          textString="*")}),
    Documentation(info="<html>
<p>
Connector with one output signal of type Real.
</p>
</html>"));
end ViewObjectNoEff2;
