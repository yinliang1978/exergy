within Exergy.Utilities;
record EnergyRateData

    Modelica.SIunits.EnergyFlowRate Ex_flow
    "Exergy flow rate from the heat flow";
    Modelica.SIunits.EnergyFlowRate E_flow;

 // Modelica.SIunits.InternalEnergy Ex_int(start=0)
   // "Exergy flow rate from the heat flow";
  //  Modelica.SIunits.InternalEnergy E_int( start=0)
  //  "Heat flow rate (positive if flowing from outside into the component)";

  //    Modelica.SIunits.EntropyFlowRate S_flow=Q/T
  //    "Entropy flow rate from the heat flow";
//equation
 //  Ex=der(Ex_int);
  // E=der(E_int);
end EnergyRateData;
