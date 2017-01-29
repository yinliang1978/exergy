within Exergy.XBuildings.Media.Water;
block specificExergyBlock "Specific Exergy"

 input ThermodynamicState state "Thermodynamic state record";
 input Exergy.Utilities.RefEnv refEnv "Ref state";
 input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.Flow;
 input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
 output Exergy.XBuildings.Media.SpecificExergy ex "Specific entropy";
// output Real dh;
// output Real ds;
// output Real dchem;
// Real h,cp, Y[2], refY[2], s;
 /*
    Modelica.Media.Interfaces.Types.SaturationProperties sata=Modelica.Media.Interfaces.Types.SaturationProperties(
                             Exergy.XModelica.Media.Water.StandardWater.saturationPressure(refEnv.T), refEnv.T);
   ThermodynamicState sat(p=sata.psat,T=sata.Tsat);
   ThermodynamicState middleState=sat;
   ThermodynamicState finalState=setState_pTX(refEnv.p,refEnv.T);

   Real h,h_c,s,s_c,chem;
   Real flowEx,volEx;
   */

Real u,s,chem;
Real flowEx,volEx;

equation
 // ex:=specificEnthalpy(state);
 u = cv_const*(state.T - refEnv.T);
 s = cv_const*Modelica.Math.log(state.T/refEnv.T);
 //  chem:=-Modelica.Constants.R*refEnv.T*(1/fluidConstants[1].molarMass)*Modelica.Math.log(max(refEnv.phi, Modelica.Constants.eps));
if eqmType==Exergy.Utilities.Types.EqmTypes.TMC then
chem= -Modelica.Constants.R*refEnv.T*(1/fluidConstants[1].molarMass)*
     Modelica.Math.log((refEnv.phi));
else
  chem= 0;
end if;
volEx= u - refEnv.T*s + chem;
// volEx:=flowEx+(refEnv.p-state.p)/(density(state));
flowEx= volEx - (refEnv.p - state.p)/d_const;

 if exType == Exergy.Utilities.Types.ExTypes.Flow then
   ex= flowEx;
 elseif exType == Exergy.Utilities.Types.ExTypes.Vol then
   ex= volEx;
 else
    ex= 0;
         Modelica.Utilities.Streams.error(
         "The exergy type " + String(exType) +
         "does not exist ");
 end if;

end specificExergyBlock;
