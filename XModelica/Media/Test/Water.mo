within Exergy.XModelica.Media.Test;
model Water
  //Exergy p105
 // Modelica.SIunits.Temperature T=273.15+2.2;
 // Modelica.SIunits.Pressure p=98.80*1000;
 // Real specificHum=0.002;
 // Real molSpecificHum=1.608*specificHum;
 // Modelica.SIunits.MoleFraction vaporMolFra=molSpecificHum/(1+molSpecificHum);
   // Exergy.Utilities.RefEnv refEnv(p=p,T=T,vaporMolFra=specificHumToVaporMolFra(specificHum));
  Exergy.Utilities.RefEnv refEnv;
Exergy.XModelica.Media.Water.StandardWater.BaseProperties bro;
Exergy.XBuildings.Media.SpecificExergy ex_Flow_TMC;
Exergy.XBuildings.Media.SpecificExergy ex_Vol_TMC;
Exergy.XBuildings.Media.SpecificExergy ex_Flow_TM;
Exergy.XBuildings.Media.SpecificExergy ex_Vol_TM;
    output Real dh;
    output Real ds;
    output Real dchem;
  //  Exergy.Utilities.RefEnv refEnv;

  // Real cp;
 // Real[2] Y=Exergy.XBuildings.Media.Air.moleToMassFractions({0.01,1 -
 //    0.01}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
 //    Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});
  //  Real[2] Y=Exergy.XBuildings.Media.Air.moleToMassFractions({0.00672,1 -
  //      0.00672}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
//    Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});

equation
bro.state.p =  98.80*1000+100000*time;
bro.state.T = 273.15+22;
//bro.state.p = 500*1000;
//bro.state.T = 273.15+13.1;
//ex=18541.3;
(ex_Flow_TMC,dh,ds,dchem) = Exergy.XModelica.Media.Water.StandardWater.specificExergy(bro.state,
  refEnv,
  Exergy.Utilities.Types.ExTypes.Flow,
  Exergy.Utilities.Types.EqmTypes.TMC);
  ex_Vol_TMC  = Exergy.XModelica.Media.Water.StandardWater.specificExergy(bro.state,
  refEnv,
  Exergy.Utilities.Types.ExTypes.Vol,
  Exergy.Utilities.Types.EqmTypes.TMC);

    ex_Flow_TM  = Exergy.XModelica.Media.Water.StandardWater.specificExergy(bro.state,
  refEnv,
  Exergy.Utilities.Types.ExTypes.Flow,
  Exergy.Utilities.Types.EqmTypes.TM);
    ex_Vol_TM  = Exergy.XModelica.Media.Water.StandardWater.specificExergy(bro.state,
  refEnv,
  Exergy.Utilities.Types.ExTypes.Vol,
  Exergy.Utilities.Types.EqmTypes.TM);

  //Exergy.Utilities.Types.ExTypes.InternalEnergyExergy
  // cp = Exergy.XBuildings.Media.Air.specificHeatCapacityCp(bro.state);

end Water;
