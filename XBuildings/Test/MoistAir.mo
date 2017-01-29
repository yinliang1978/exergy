within Exergy.XBuildings.Test;
model MoistAir
  Exergy.XBuildings.Media.Air.BaseProperties bro;
  Exergy.XBuildings.Media.SpecificExergy ex;
  Real dex;
  Exergy.Utilities.RefEnv refEnv(T=273.15 + 15,state(X=Y));
  Real cp;

  Real[2] Y=Exergy.XBuildings.Media.Air.moleToMassFractions({0.00672,1 -
      0.00672}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
      Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});
  Real[2] Y2=Exergy.XBuildings.Media.Air.moleToMassFractions({0.0983,1 -
      0.0983}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
      Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});
equation
  bro.state.p = 101325;
  bro.state.T = 273.15 + 60;
  bro.state.X[1] = Y2[1];
  //ex=18541.3;
  ex = Exergy.XBuildings.Media.Air.specificEnthalpyExergy(bro.state,
    refEnv.state,Exergy.Utilities.Types.ExTypes.InternalEnergyExergy);
  cp = Exergy.XBuildings.Media.Air.specificHeatCapacityCp(bro.state);
  dex=der(ex);
end MoistAir;
