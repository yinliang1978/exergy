within Exergy.Contours;
function MoistPrint
  input Exergy.Utilities.RefEnv refEnv;
  input Real T;
  input Real X;
  input String fileName;
protected
  Exergy.XBuildings.Media.Air.ThermodynamicState state;
  Exergy.XBuildings.Media.Air.ThermodynamicState stateCom;
  Real ex,phi;
  Real ps,omega,k_mair,satX;
  String outs;
algorithm

  k_mair:=Modelica.Media.Air.MoistAir.k_mair;
  state.p:=refEnv.p;

 // for T in 273:1:(273+60) loop
 //   for X in 0.0001:0.001:0.042 loop
  state.T:=T;
  state.X[1]:=X;
  state.X[2]:=1-state.X[1];

  stateCom:=state;

  ps:=Modelica.Media.Air.MoistAir.saturationPressure(T);
  omega:=k_mair*ps/(state.p-ps);
  satX:=omega/(1.0+omega);

  if X>=satX then
    stateCom.X[1]:=satX;
    stateCom.X[2]:=1.0-stateCom.X[1];
  end if;
    ex:=Exergy.XBuildings.Media.Air.specificExergy(stateCom, refEnv);
    phi:=Modelica.Media.Air.MoistAir.relativeHumidity(stateCom);
  // ex:=1;
  //phi:=Exergy.XModelica.Media.Air.MoistAir.relativeHumidity_pTX(state.p,state.T,state.X);
  //if phi<0.999 then
 // ex:=Exergy.XBuildings.Media.Air.specificExergy(state, refEnv);
 //outs:=String(T)+","+String(state.X[1])+","+String(ex);
  Modelica.Utilities.Streams.print(String(T-273.15)+","+String(state.X[1]*1000/state.X[2])+","+String(ex)+","+String(phi),fileName);
//  Modelica.Utilities.Streams.close("d.txt");

  //else
  //end if;
  //  end for;
  //end for;
end MoistPrint;
