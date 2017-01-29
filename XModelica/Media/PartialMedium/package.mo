within Exergy.XModelica.Media;
partial package PartialMedium 
extends Modelica.Media.Interfaces.PartialMedium;


  replaceable function specificExergy "Specific Exergy"
    input ThermodynamicState state "Thermodynamic state record";
    input Exergy.Utilities.RefEnv refEnv "Ref state";
    input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.Flow;
    input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
    output Exergy.XBuildings.Media.SpecificExergy ex "Specific entropy";
   //  output Real dh;
     //output Real ds;
     //output Real dchem;

   //  Real ex2;
  algorithm

    ex := 0;

    annotation (derivative(noDerivative=refEnv)=specificExergy_der,Inline=true);
  end specificExergy;


      replaceable function specificExergy_der "Specific Exergy"
  input ThermodynamicState state "Thermodynamic state record";
  input Exergy.Utilities.RefEnv refEnv "Ref state";
  input Exergy.Utilities.Types.ExTypes exType=Exergy.Utilities.Types.ExTypes.Flow;
  input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
  input ThermodynamicState state_der;
  // input Real T_der;
  // input Real X_der[nX];
  //  input Real X2_der;
  output Real ex_der "Specific entropy";
  // output Real dh;
  // output Real ds;
  // output Real dchem;

  //  Real ex2;

      algorithm
      ex_der:=0;
  //   ex_der:=state_der.T;
  //ex_der := 0.89;
  // ex_der:=(1-refEnv.T/state.T)*(dryair.cp*(1 - state.X[Water]) + steam.cp*state.X[Water]);

      end specificExergy_der;

/*
 replaceable function specificExergy "Specific Exergy"
   input ThermodynamicState state "Thermodynamic state record";
   input Exergy.Utilities.RefEnv refEnv "Ref state";
    input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.Flow;
    input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
  output Exergy.XBuildings.Media.SpecificExergy ex "Specific entropy";

 algorithm 
  ex:=0;
    annotation (derivative(noDerivative=refEnv)=specificExergy_der,Inline=true);
        annotation (derivative(noDerivative=refEnv)=specificExergy_der,Inline=true);
 end specificExergy;

      replaceable function specificExergy_der "Specific Exergy"
    input ThermodynamicState state "Thermodynamic state record";
    input Exergy.Utilities.RefEnv refEnv "Ref state";
    input Exergy.Utilities.Types.ExTypes exType=Exergy.Utilities.Types.ExTypes.Flow;
    input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
    input ThermodynamicState state_der;
    // input Real T_der;
    // input Real X_der[nX];
    //  input Real X2_der;
    output Real ex_der "Specific entropy";
    // output Real dh;
    // output Real ds;
    // output Real dchem;
      algorithm 
    ex_der := 0.0;
    //   ex_der:=state_der.T;
    //ex_der := 0.89;
    // ex_der:=(1-refEnv.T/state.T)*(dryair.cp*(1 - state.X[Water]) + steam.cp*state.X[Water]);

end specificExergy_der;
*/

 /*
 replaceable function specificEnthalpyExergy_der "Specific Exergy"
  input ThermodynamicState state "Thermodynamic state record";
  input ThermodynamicState refState "Ref state";
  input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.EnthalpyExergy;
  input PotentialEnergy E=0 
      " energy from potential ,e.g. kinetic,potential due to external fields ";
  input Real dE(unit="J/(kg.s)");
  output Real ex_der(unit="J/(kg.s)") "Specific entropy";
  annotation (derivative=specificEnthalpyExergy_der,Inline=true);
  
  end specificEnthalpyExergy_der;
  */
end PartialMedium;
