within Exergy.XBuildings.Media;
package Air 
  extends Buildings.Media.Air;

  extends Modelica.Icons.Package;


  replaceable function specificExergy "Specific Exergy"
    input ThermodynamicState state "Thermodynamic state record";
    input Exergy.Utilities.RefEnv refEnv "Ref state";
    input Exergy.Utilities.Types.ExTypes exType = Exergy.Utilities.Types.ExTypes.Flow;
    input Exergy.Utilities.Types.EqmTypes eqmType=Exergy.Utilities.Types.EqmTypes.TMC;
    output Exergy.XBuildings.Media.SpecificExergy ex "Specific entropy";
   //  output Real dh;
     //output Real ds;
     //output Real dchem;
protected
    Real h,cp, Y[2], refY[2], s;
     Real omega;
   //  Real ex2;
  algorithm
      omega:=state.X[Water]/(1 - state.X[Water]);
      /*
    ex2:=(dryair.cp+omega*steam.cp)*refEnv.T*(state.T/refEnv.T-1-Modelica.Math.log(state.T/refEnv.T))
         +(1+1.6078*omega)*Modelica.Constants.R/dryair.MM*refEnv.T*Modelica.Math.log(state.p/refEnv.p)
         +Modelica.Constants.R/dryair.MM*refEnv.T*((1+1.6078*omega)*Modelica.Math.log((1+1.6078*refEnv.omega)/(1+1.6078*omega))
     +1.6078*omega*Modelica.Math.log(omega/refEnv.omega));
*/
   // h := (state.T - refState.T)*dryair.cp*(1 - state.X[Water]) + ((state.T
   //    - refState.T)*steam.cp)*state.X[Water];
    //cv := dryair.cv*(1 - state.X[Water]) + steam.cv*state.X[Water];
    cp := dryair.cp*(1 - state.X[Water]) + steam.cp*state.X[Water];

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

     h:=cp*(state.T - refEnv.T);

      Y := massToMoleFractions(state.X, {steam.MM,dryair.MM});
      //refY := refEnv.MoistMolX;
      refY :=massToMoleFractions(refEnv.MoistMassX, {steam.MM,dryair.MM});

    s := cp*Modelica.Math.log(state.T/refEnv.T) - (Modelica.Constants.R
      /molarMass(state))*Modelica.Math.log(state.p/refEnv.p);

    ex := h - refEnv.T*s + (Modelica.Constants.R/molarMass(state))*
    refEnv.T*sum(Y[i]*Modelica.Math.log(max(Y[i]/refY[i], Modelica.Constants.eps))
      for i in 1:2);
            //ex := h - refState.T*s;
  if exType == Exergy.Utilities.Types.ExTypes.Flow then
    ex:=ex;
  elseif exType == Exergy.Utilities.Types.ExTypes.Vol then
    ex:=ex+(refEnv.p-state.p)/(density(state));
  else
          Modelica.Utilities.Streams.error(
          "The exergy type " + String(exType) +
          "does not exist ");
  end if;

  //    dh:=h;
   // ds:=s;
  //  dchem:=(1.0/Exergy.XModelica.Media.Air.MoistAir.k_mair);
  //  dchem:=Modelica.Math.log(max(Y[1]/refY[1], Modelica.Constants.eps));
    //   dchem:=(Modelica.Constants.R/molarMass(state))*
    //   refEnv.T*sum(Y[i]*Modelica.Math.log(max(Y[i]/refY[i], Modelica.Constants.eps))
  //  for i in 1:2);

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
  Real Y[2], refY[2], MM_a, MM_v, y_v, y_0v, C_pa, C_pv, T, T_0, R, P,
      P_0,
      k_a;
  Real x_v, C_p;
  Real MM, ex;
  Real exDy_v, exDT, y_vDx_v;
  //  Real ex2;

      algorithm
      Y := massToMoleFractions(state.X, {steam.MM,dryair.MM});
      //refY := refEnv.MoistMolX;
      refY :=massToMoleFractions(refEnv.MoistMassX, {steam.MM,dryair.MM});

      MM_a:=dryair.MM;
      MM_v:=steam.MM;
      y_v:=Y[1];
      y_0v:=refY[1];
      C_pa:=dryair.cp*MM_a;
      C_pv:=steam.cp*MM_v;
      T:=state.T;
      T_0:=refEnv.T;
      R:=Modelica.Constants.R;
      P:=state.p;
      P_0:=refEnv.p;
      k_a:=MM_a/MM_v;

      MM:=MM_a*(1-y_v)+MM_v*y_v;
      //C_p:=(C_pa*(1-state.X[Water])+C_pv*state.X[Water])*MM;
      C_p:=C_pa*(1-y_v)+C_pv*y_v;
      ex:=R*T_0*(y_v*Modelica.Math.log(y_v/y_0v)+(1-y_v)*Modelica.Math.log((1-y_v)/(1-y_0v)))
         +R*T_0*Modelica.Math.log(P/P_0)
         +C_p*(T-T_0-T_0*Modelica.Math.log(T/T_0));

      exDy_v:=(MM_a-MM_v)*ex/(MM*MM)
              +(R*T_0*(Modelica.Math.log(y_v/y_0v)-Modelica.Math.log((1-y_v)/(1-y_0v)))
                +(C_pv-C_pa)*(T-T_0-T_0*Modelica.Math.log(T/T_0)))/MM;

      exDT:=(1-T_0/T)*C_p/MM;

      x_v:=state.X[1];
      y_vDx_v:=k_a/((1+(k_a-1)*x_v)*(1+(k_a-1)*x_v));

      ex_der:=exDy_v*y_vDx_v*state_der.X[1]+exDT*state_der.T;
  //   ex_der:=state_der.T;
  //ex_der := 0.89;
  // ex_der:=(1-refEnv.T/state.T)*(dryair.cp*(1 - state.X[Water]) + steam.cp*state.X[Water]);

      end specificExergy_der;


end Air;
