within Exergy.XModelica.Media.Test;
model MoistAir

    //Exergy p105
 // Modelica.SIunits.Temperature T=273.15+2.2;
 // Modelica.SIunits.Pressure p=98.80*1000;
 // Real specificHum=0.002;
 // Real molSpecificHum=1.608*specificHum;
 // Modelica.SIunits.MoleFraction vaporMolFra=molSpecificHum/(1+molSpecificHum);
//   Exergy.Utilities.RefEnv refEnv(p=p,T=T,vaporMolFra=specificHumToVaporMolFra(specificHum));

 //  Exergy.Utilities.RefEnv refEnv;

   Modelica.SIunits.Temperature T=273.15+2.2;
 // Modelica.SIunits.Pressure p=98.80*1000;
   Real specificHum=0.002;
    Exergy.Utilities.RefEnv refEnv(T=T,vaporMolFra=specificHumToVaporMolFra(specificHum));

   Exergy.XModelica.Media.Air.MoistAir.BaseProperties bro;

     Real phi=Exergy.XModelica.Media.Air.MoistAir.relativeHumidity_pTX(bro.state.p,bro.state.T, bro.state.X);
Exergy.XBuildings.Media.SpecificExergy ex;
//Real dex;
//Exergy.Utilities.RefEnv refEnv(T=273.15 + 15,state(X=Y));
//Exergy.Utilities.RefEnv refEnv;
//Real cp;

//Real[2] Y=Exergy.XBuildings.Media.Air.moleToMassFractions({0.00672,1 -
//    0.00672}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
//    Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});
//Real[2] Y2=Exergy.XBuildings.Media.Air.moleToMassFractions({0.0183,1 -
//    0.0183}, {Modelica.Media.IdealGases.Common.SingleGasesData.H2O.MM,
//    Modelica.Media.IdealGases.Common.SingleGasesData.Air.MM});
Real[2] Y2={0.001,1.0-0.01};
    output Real dh;
    output Real ds;
    output Real dchem;
equation
bro.state.p = 101325;
bro.state.T = 273.15 + 14;
bro.state.X[1] = Exergy.XModelica.Media.Test.specificHumToVaporMassFra(0.004);
//ex=18541.3;
(ex,dh,ds,dchem) = Exergy.XModelica.Media.Air.MoistAir.specificExergy(bro.state,
  refEnv,Exergy.Utilities.Types.ExTypes.Flow);
//cp = Exergy.XBuildings.Media.Air.specificHeatCapacityCp(bro.state);
//dex=der(ex);
end MoistAir;
