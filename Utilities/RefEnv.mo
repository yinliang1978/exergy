within Exergy.Utilities;
record RefEnv
  extends Modelica.Icons.Record;

  //replaceable package AirMedium = Buildings.Media.Air "Medium in the component";
  replaceable package MoistMedium =Modelica.Media.Air.MoistAir;
  //temperature and pressure
  Modelica.SIunits.Temperature T=T0;
  Modelica.SIunits.Pressure p=p0;

  //must from vaperMolFra setting ,else refValue is error

  Modelica.SIunits.MoleFraction vaporMolFra=airMolFra[AirComposition.H2O];
  Modelica.SIunits.MoleFraction dryAirMolFra=1.0-vaporMolFra;
  Modelica.SIunits.MoleFraction MoistMolX[2]={vaporMolFra,dryAirMolFra};
  Modelica.SIunits.MolarMass MoistMM=MoistMolX*MoistMedium.MMX;

  Modelica.SIunits.MassFraction MoistMassX[2]=
                   MoistMedium.moleToMassFractions(MoistMolX,MoistMedium.MMX);
  Real omega=MoistMassX[1]/MoistMassX[2];

  Real phi=MoistMedium.relativeHumidity_pTX(p,T, MoistMassX);
    /*
  Modelica.SIunits.MassFraction vaporMassFra=airMolFra[AirComposition.H2O];
  Modelica.SIunits.MassFraction dryAirMassFra=1.0-vaporMassFra;
  Modelica.SIunits.MassFraction MoistX[2]={vaporMassFra,dryAirMassFra};
  Modelica.SIunits.MolarMass MoistMM=
     1/(MoistX[MoistMedium.Water]/MoistMedium.MMX[MoistMedium.Water]
     + (MoistX[MoistMedium.Air])/MoistMedium.MMX[MoistMedium.Air]);
//  MM = 1/(Xi[Water]/MMX[Water] + (1.0 - Xi[Water])/MMX[Air]);
Modelica.SIunits.MoleFraction vaporMolFra=airMolFra[AirComposition.H2O];

*/
  //moleToMassFractions
  MoistMedium.ThermodynamicState state(
    p=p,
    T=T,
    X=MoistMassX);

  constant Modelica.SIunits.Temperature T0=298.15;

 // constant Modelica.SIunits.Pressure p0=101325;
  constant Modelica.SIunits.Pressure p0=1e5;

  type AirComposition = enumeration(
      N2,
      O2,
      H2O,
      Ar,
      CO2,
      Ne,
      He) "Enumeration defining the reference saturated air" annotation (
      Evaluate=true);

       constant Modelica.SIunits.MoleFraction airMolFra[AirComposition]={0.7561,
  0.2028,0.03167,0.0091,0.0003,1.77e-5,5.08e-6};

      /*
  constant Modelica.SIunits.MoleFraction airMolFra[AirComposition]={0.7567,
  0.2035,0.0303,0.0091,0.0003,0.0001};
  */

  type ConComposition = enumeration(
      H2O,
      CaCO3,
      CaSO4D2H2O) "Enumeration defining the condensed phases" annotation (
      Evaluate=true);
/*
  type ReferenceStateOfWater = enumeration(
      LiquidWater,
      MoistAir) 
    "Enumeration defining Reference State Of Water when evaluating exergy"
    annotation (Evaluate=true);
    */

end RefEnv;
