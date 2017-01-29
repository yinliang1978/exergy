within Exergy.XBuildings.Media.Air;
block specificExergyBlock

    input ThermodynamicState state "Thermodynamic state record";
input Exergy.Utilities.RefEnv refEnv "Ref state";
input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.Flow;
input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
output Exergy.XBuildings.Media.SpecificExergy ex "Specific entropy";
// output Real dh;
// output Real ds;
// output Real dchem;
//protected
Real h,cp, Y[2], refY[2], s;
 Real omega,extemp;
  //  Real ex2;
equation
  omega=state.X[Water]/(1 - state.X[Water]);
  /*
    ex2:=(dryair.cp+omega*steam.cp)*refEnv.T*(state.T/refEnv.T-1-Modelica.Math.log(state.T/refEnv.T))
         +(1+1.6078*omega)*Modelica.Constants.R/dryair.MM*refEnv.T*Modelica.Math.log(state.p/refEnv.p)
         +Modelica.Constants.R/dryair.MM*refEnv.T*((1+1.6078*omega)*Modelica.Math.log((1+1.6078*refEnv.omega)/(1+1.6078*omega))
     +1.6078*omega*Modelica.Math.log(omega/refEnv.omega));
*/
  // h := (state.T - refState.T)*dryair.cp*(1 - state.X[Water]) + ((state.T
  //    - refState.T)*steam.cp)*state.X[Water];
//cv := dryair.cv*(1 - state.X[Water]) + steam.cv*state.X[Water];
cp = dryair.cp*(1 - state.X[Water]) + steam.cp*state.X[Water];

  //  if exType == Exergy.Utilities.Types.ExTypes.EnthalpyExergy then
  //    h:=cp*(state.T - refState.T);
  //  elseif exType == Exergy.Utilities.Types.ExTypes.InternalEnergyExergy then
 // h:=cp*(state.T - refState.T)- (pStp/dStp-pStp/dStp)+refState.p*(1.0/(state.p*dStp/pStp)-1.0/(refState.p*dStp/pStp));
  //   h:=cp*(state.T - refState.T)+pStp/dStp*(refState.p/(state.p)-1.0);

//   else
  //   Modelica.Utilities.Streams.error(
 //     "The exergy type " + String(exType) +
   //   "does not exist ");
  // end if;

 h=cp*(state.T - refEnv.T);

  Y = massToMoleFractions(state.X, {steam.MM,dryair.MM});
  //refY := refEnv.MoistMolX;
  refY =massToMoleFractions(refEnv.MoistMassX, {steam.MM,dryair.MM});

s = cp*Modelica.Math.log(state.T/refEnv.T) - (Modelica.Constants.R
  /molarMass(state))*Modelica.Math.log(state.p/refEnv.p);

extemp = h - refEnv.T*s + (Modelica.Constants.R/molarMass(state))*
refEnv.T*sum(Y[i]*Modelica.Math.log(max(Y[i]/refY[i], Modelica.Constants.eps))
  for i in 1:2);
        //ex := h - refState.T*s;
  if exType == Exergy.Utilities.Types.ExTypes.Flow then
ex=extemp;
  elseif exType == Exergy.Utilities.Types.ExTypes.Vol then
ex=extemp+(refEnv.p-state.p)/(density(state));
  else
    ex=0;
      Modelica.Utilities.Streams.error(
      "The exergy type " + String(exType) +
      "does not exist ");
  end if;

 // dh:=h;
//ds:=s;
  //  dchem:=(1.0/Exergy.XModelica.Media.Air.MoistAir.k_mair);
  //  dchem:=Modelica.Math.log(max(Y[1]/refY[1], Modelica.Constants.eps));
//dchem:=(Modelica.Constants.R/molarMass(state))*
//refEnv.T*sum(Y[i]*Modelica.Math.log(max(Y[i]/refY[i], Modelica.Constants.eps))
//for i in 1:2);

end specificExergyBlock;
