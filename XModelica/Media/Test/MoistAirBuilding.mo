within Exergy.XModelica.Media.Test;
model MoistAirBuilding

    //Exergy p105
  // Modelica.SIunits.Temperature T=273.15+0;
 // Modelica.SIunits.Pressure p=98.80*1000;
  // Real specificHum=0.006425;
 // Real molSpecificHum=1.608*specificHum;
 // Modelica.SIunits.MoleFraction vaporMolFra=molSpecificHum/(1+molSpecificHum);
 //   Exergy.Utilities.RefEnv refEnv(p=p,T=T,vaporMolFra=specificHumToVaporMolFra(specificHum));

  Real omega=0.000718;
  Real MoistMassX[2]={omega/(1+omega),1/(1+omega)};
  Exergy.Utilities.RefEnv refEnv(T=T,MoistMassX=MoistMassX);

    Modelica.SIunits.Temperature T=273.15-9.9;
// Modelica.SIunits.Pressure p=98.80*1000;
   Real specificHum=0.01+0.02*time;
  //    Exergy.Utilities.RefEnv refEnv(T=T,vaporMolFra=specificHumToVaporMolFra(specificHum));

    Exergy.XBuildings.Media.Air.BaseProperties bro;

     Real phi=Exergy.XModelica.Media.Air.MoistAir.relativeHumidity_pTX(bro.state.p,bro.state.T, bro.state.X);
     Exergy.XBuildings.Media.SpecificExergy ex;
    // input  Real ex2;
    Real ex2;
     Exergy.XBuildings.Media.Air.specificExergyBlock exB(state=(bro.state),refEnv=refEnv);
      Real dex,dex2;
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
  //   output Real dh;
  //  output Real ds;
  //  output Real dchem;
equation
 // connect(ex2,exB.ex);
 ex2=exB.ex;
bro.state.p = 100000;
bro.state.T = 273.15 +60+time;
bro.state.X[1] = Exergy.XModelica.Media.Test.specificHumToVaporMassFra(specificHum);
//ex=18541.3;
/*
(ex,dh,ds,dchem) =  Exergy.XBuildings.Media.Air.specificExergy(bro.state,
refEnv,Exergy.Utilities.Types.ExTypes.Flow);
*/
ex  =  Exergy.XBuildings.Media.Air.specificExergy(bro.state,
refEnv,Exergy.Utilities.Types.ExTypes.Flow);
dex=der(ex);
dex2=der(ex2);
//cp = Exergy.XBuildings.Media.Air.specificHeatCapacityCp(bro.state);
//dex=der(ex);
end MoistAirBuilding;
