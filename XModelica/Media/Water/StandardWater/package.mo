within Exergy.XModelica.Media.Water;
package StandardWater 
  extends Modelica.Media.Water.StandardWater;

  replaceable function specificExergy "Specific Exergy"
      input ThermodynamicState state "Thermodynamic state record";
      input Exergy.Utilities.RefEnv refEnv "Ref state";
      input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.Flow;
      input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
      output Exergy.XBuildings.Media.SpecificExergy ex "Specific entropy";
      output Real dh;
      output Real ds;
      output Real dchem;
     // Real h,cp, Y[2], refY[2], s;
protected
     ThermodynamicState sat=setBubbleState(SaturationProperties(saturationPressure(refEnv.T), refEnv.T));
     ThermodynamicState middleState=sat;
     ThermodynamicState finalState=setState_pTX(refEnv.p,refEnv.T);

     Real h,h_c,s,s_c,chem;
     Real flowEx,volEx;

  algorithm
    h:=specificEnthalpy(state);
    s:=specificEntropy(state);
    if eqmType==Exergy.Utilities.Types.EqmTypes.TMC then
    h_c:=specificEnthalpy(middleState);
    s_c:=specificEntropy(middleState);
    chem:=-Modelica.Constants.R*refEnv.T*(1/fluidConstants[1].molarMass)*Modelica.Math.log(max(refEnv.phi, Modelica.Constants.eps));
    flowEx:=h-h_c-refEnv.T*(s-s_c)+(state.p-middleState.p)/(density(middleState))+chem;
    else
    h_c:=specificEnthalpy(finalState);
    s_c:=specificEntropy(finalState);
    chem:=0;
    flowEx:=h-h_c-refEnv.T*(s-s_c);
    end if;

    //chem:=0;
     // for energy is only the function of temperature,the following flowEx is as volEx
     //consider the
    //volEx:=h-h_c-refEnv.T*(s-s_c)+chem;
    // flowEx:=volEx-(refEnv.p-state.p)/(density(state));

     // flowEx:=h-h_c-refEnv.T*(s-s_c)+(state.p-middleState.p)/(density(middleState))-chem;
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
  // dchem:=specificHeatCapacityCp(state);

  end specificExergy;
end StandardWater;
