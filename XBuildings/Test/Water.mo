within Exergy.XBuildings.Test;
model Water
  Exergy.XBuildings.Media.Water.BaseProperties bro;
  Exergy.XBuildings.Media.SpecificExergy ex;
 //  Exergy.Utilities.RefEnv refEnv;
   Exergy.Utilities.RefEnv refEnv(state(X=Y));
 // Real cp;
    Real[2] Y=Exergy.XBuildings.Media.Air.moleToMassFractions({0.01,1 -
       0.01}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
       Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});
//  Real[2] Y=Exergy.XBuildings.Media.Air.moleToMassFractions({0.00672,1 -
//      0.00672}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
  //    Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});

equation
  bro.state.p = 101325;
  bro.state.T = 298.15;
  //ex=18541.3;
  ex = Exergy.XBuildings.Media.Water.specificEnthalpyExergy(bro.state,
    refEnv.state,
    Exergy.Utilities.Types.ExTypes.EnthalpyExergy);
    //Exergy.Utilities.Types.ExTypes.InternalEnergyExergy
 // cp = Exergy.XBuildings.Media.Air.specificHeatCapacityCp(bro.state);

end Water;
