within Exergy.XModelica.Media.Air;
package MoistAir 
  extends Modelica.Media.Air.MoistAir;

  replaceable function specificExergy "Specific Exergy"
    input ThermodynamicState state "Thermodynamic state record";
    input Exergy.Utilities.RefEnv refEnv "Ref state";
    input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.Flow;
    input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
    output Exergy.XBuildings.Media.SpecificExergy ex "Specific entropy";
   // Real h,cp, Y[2], refY[2], s;

       output Real dh;
      output Real ds;
      output Real dchem;
protected
   ThermodynamicState middleState(p=refEnv.p,T=refEnv.T,X=state.X);
   ThermodynamicState finalState(p=refEnv.p,T=refEnv.T,X=refEnv.MoistMassX);
   Modelica.SIunits.MoleFraction MoistMolX[2]=massToMoleFractions(state.X,MMX);
   Modelica.SIunits.MolarMass MoistMM=MoistMolX*MMX;
   Real h,h_c,s,s_c,chem;
   Real flowEx,volEx;

  algorithm
  h:=Exergy.XModelica.Media.Air.MoistAir.specificEnthalpy(state);
  h_c:=Exergy.XModelica.Media.Air.MoistAir.specificEnthalpy(middleState);
  s:=Exergy.XModelica.Media.Air.MoistAir.specificEntropy(state);
  s_c:=Exergy.XModelica.Media.Air.MoistAir.specificEntropy(middleState);

  chem:=Modelica.Constants.R*refEnv.T*(1/MoistMM)*
        sum({MoistMolX[i]*Modelica.Math.log(max(MoistMolX[i]/refEnv.MoistMolX[i], Modelica.Constants.eps))
        for i in 1:2});
  //chem:=0;

    flowEx:=h-h_c-refEnv.T*(s-s_c)+chem;
  //  flowEx := h - refEnv.T*s;
  volEx:=flowEx+(refEnv.p-state.p)/(density(state));

  if exType == Exergy.Utilities.Types.ExTypes.Flow then
    ex:=flowEx;
  elseif exType == Exergy.Utilities.Types.ExTypes.Vol then
    ex:=volEx;
  else
          Modelica.Utilities.Streams.error(
          "The exergy type " + String(exType) +
          "does not exist ");
  end if;
      dh:=h-h_c;
      ds:=s-s_c;
      dchem:=chem;

  /*
   // h := (state.T - refState.T)*dryair.cp*(1 - state.X[Water]) + ((state.T
   //    - refState.T)*steam.cp)*state.X[Water];
    //cv := dryair.cv*(1 - state.X[Water]) + steam.cv*state.X[Water];
    cp := dryair.cp*(1 - state.X[Water]) + steam.cp*state.X[Water];

    if exType == Exergy.Utilities.Types.ExTypes.EnthalpyExergy then
      h:=cp*(state.T - refState.T);
    elseif exType == Exergy.Utilities.Types.ExTypes.InternalEnergyExergy then
     // h:=cp*(state.T - refState.T)- (pStp/dStp-pStp/dStp)+refState.p*(1.0/(state.p*dStp/pStp)-1.0/(refState.p*dStp/pStp));
      h:=cp*(state.T - refState.T)+pStp/dStp*(refState.p/(state.p)-1.0);

    else
      Modelica.Utilities.Streams.error(
          "The exergy type " + String(exType) +
          "does not exist ");
    end if;

    Y := massToMoleFractions(state.X, {steam.MM,dryair.MM});
    refY := massToMoleFractions(refState.X, {steam.MM,dryair.MM});

    s := cp*Modelica.Math.log(state.T/refState.T) - (Modelica.Constants.R
      /molarMass(state))*Modelica.Math.log(state.p/refState.p);

    ex := h - refState.T*s + (Modelica.Constants.R/molarMass(state))*
    refState.T*sum(Y[i]*Modelica.Math.log(max(Y[i]/refY[i], Modelica.Constants.eps))
      for i in 1:2);
      //ex := h - refState.T*s;
      
      */

   // annotation (derivative=specificEnthalpyExergy_der,Inline=true);
  end specificExergy;
end MoistAir;
