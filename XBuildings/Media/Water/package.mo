within Exergy.XBuildings.Media;
package Water 
  extends Buildings.Media.Water;
  // the density is constant


  replaceable function specificExergy "Specific Exergy"

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
protected
   Real u,s,chem;
   Real flowEx,volEx;

  algorithm
    // ex:=specificEnthalpy(state);
    u := cv_const*(state.T - refEnv.T);
    s := cv_const*Modelica.Math.log(state.T/refEnv.T);
  //  chem:=-Modelica.Constants.R*refEnv.T*(1/fluidConstants[1].molarMass)*Modelica.Math.log(max(refEnv.phi, Modelica.Constants.eps));
   if eqmType==Exergy.Utilities.Types.EqmTypes.TMC then
   chem:=-Modelica.Constants.R*refEnv.T*(1/fluidConstants[1].molarMass)*Modelica.Math.log((refEnv.phi));
   else
     chem:=0;
   end if;
   volEx:= u- refEnv.T*s+chem;
   // volEx:=flowEx+(refEnv.p-state.p)/(density(state));
   flowEx:=volEx-(refEnv.p-state.p)/d_const;

    if exType == Exergy.Utilities.Types.ExTypes.Flow then
      ex:=flowEx;
    elseif exType == Exergy.Utilities.Types.ExTypes.Vol then
      ex:=volEx;
    else
            Modelica.Utilities.Streams.error(
            "The exergy type " + String(exType) +
            "does not exist ");
    end if;

      /*
  
  h:=specificEnthalpy(state);
  s:=specificEntropy(state);
  if eqmType==Exergy.Utilities.Types.EqmTypes.TMC then
  h_c:=specificEnthalpy(middleState);
  s_c:=specificEntropy(middleState);
  chem:=-Modelica.Constants.R*refEnv.T*(1/fluidConstants[1].molarMass)*Modelica.Math.log(max(refEnv.phi, Modelica.Constants.eps));
  else
  h_c:=specificEnthalpy(finalState);
  s_c:=specificEntropy(finalState);
  chem:=0;
  end if;

    //chem:=0;
     // for energy is only the function of temperature,the following flowEx is as volEx
     //consider the
     volEx:=h-h_c-refEnv.T*(s-s_c)+chem;
     flowEx:=volEx-(refEnv.p-state.p)/(density(state));
    //  flowEx:=h-h_c-refEnv.T*(s-s_c)+chem;
    //  flowEx := h - refEnv.T*s;
    //volEx:=flowEx+(refEnv.p-state.p)/(density(state));

    if exType == Exergy.Utilities.Types.ExTypes.Flow then
      ex:=flowEx;
    elseif exType == Exergy.Utilities.Types.ExTypes.Vol then
      ex:=volEx;
    else
            Modelica.Utilities.Streams.error(
            "The exergy type " + String(exType) +
            "does not exist ");
    end if;
    
    
    
    */
  /*
  dh:=h-h_c-(refEnv.p-state.p)/(density(state));
  ds:=s-s_c;
  dchem:=chem;
  */

  /*
    input ThermodynamicState state "Thermodynamic state record";
    input Air.ThermodynamicState refState "Ref state";
    input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.EnthalpyExergy;
    input PotentialEnergy E=0 
      " energy from potential ,e.g. kinetic,potential due to external fields ";
    output SpecificExergy ex "Specific entropy";
  protected 
    Real h,cp,cv, s;
  algorithm 
   // h := (state.T - refState.T)*dryair.cp*(1 - state.X[Water]) + ((state.T
   //    - refState.T)*steam.cp)*state.X[Water];
    //cv := dryair.cv*(1 - state.X[Water]) + steam.cv*state.X[Water];
    cp := cp_const;
    cv := cv_const;
    if exType == Exergy.Utilities.Types.ExTypes.EnthalpyExergy then
      h:=cp*(state.T - refState.T)+( state.p-Exergy.XBuildings.Media.Air.saturationPressure(state.T)) /d_const;
    elseif exType == Exergy.Utilities.Types.ExTypes.InternalEnergyExergy then
     // h:=cp*(state.T - refState.T)- (pStp/dStp-pStp/dStp)+refState.p*(1.0/(state.p*dStp/pStp)-1.0/(refState.p*dStp/pStp));
      h:=cv*(state.T - refState.T);
     elseif exType == Exergy.Utilities.Types.ExTypes.LiquidEnthalpyExergy then
       h:=cp*(state.T - refState.T)+( state.p-refState.p) /d_const;
       elseif exType == Exergy.Utilities.Types.ExTypes.LiquidInternalEnergyExergy then
         h:=cv*(state.T - refState.T);
    else
      Modelica.Utilities.Streams.error(
          "The exergy type " + String(exType) +
          "does not exist ");
    end if;

    s :=cv_const*Modelica.Math.log(state.T/refState.T);

      if exType == Exergy.Utilities.Types.ExTypes.LiquidEnthalpyExergy or exType == Exergy.Utilities.Types.ExTypes.LiquidInternalEnergyExergy then
        ex := h - refState.T*s;
      else
        ex := h - refState.T*s
         - (Modelica.Constants.R/molarMass(state))*
        refState.T*Modelica.Math.log(max(Modelica.Media.Air.MoistAir.relativeHumidity(refState),Modelica.Constants.eps));
        end if
        */
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
protected
     Real C_v,T,T_0,V_0;
     Real exDp,exDT;
    //  Real ex2;

   algorithm
     C_v:=cv_const;
     T:=state.T;
     T_0:=refEnv.T;
     V_0:=1.0/d_const;
     exDT:=C_v-T_0*C_v/T;
     exDp:=V_0;
    if exType == Exergy.Utilities.Types.ExTypes.Flow then

      ex_der:=state_der.p*exDp+state_der.T*exDT;
    elseif exType == Exergy.Utilities.Types.ExTypes.Vol then
      ex_der:=state_der.T*exDT;
    else
            Modelica.Utilities.Streams.error(
            "The exergy type " + String(exType) +
            "does not exist ");
    end if;
    //   ex_der:=state_der.T;
    //ex_der := 0.89;
    // ex_der:=(1-refEnv.T/state.T)*(dryair.cp*(1 - state.X[Water]) + steam.cp*state.X[Water]);

   end specificExergy_der;


end Water;
