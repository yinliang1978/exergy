within Exergy.Utilities;
record EnergyData

  //  Modelica.SIunits.EnergyFlowRate Ex "Exergy flow rate from the heat flow";
  //  Modelica.SIunits.EnergyFlowRate E;

   Modelica.SIunits.InternalEnergy Ex(start=0,stateSelect=StateSelect.avoid)
    "Exergy flow rate from the heat flow";
     Modelica.SIunits.InternalEnergy E( start=0,stateSelect=StateSelect.avoid)
    "Heat flow rate (positive if flowing from outside into the component)";

  //    Modelica.SIunits.EntropyFlowRate S_flow=Q/T
  //    "Entropy flow rate from the heat flow";
//equation
 //  Ex=der(Ex_int);
  // E=der(E_int);
end EnergyData;
