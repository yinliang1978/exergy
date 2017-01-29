within Exergy.Contours;
function WaterPrint_ph
  input Exergy.Utilities.RefEnv refEnv;
  input Real p;
  input Real h;
  input String fileName;

protected
  Exergy.XModelica.Media.Water.StandardWater.ThermodynamicState state;
 // Exergy.XModelica.Media.Water.StandardWater.ThermodynamicState stateCom;
   Real ex,T,vaperQ;//,phi;
 // Real ps,omega,k_mair,satX;
//  String outs;
algorithm
  //setState_dTX
  state:=Exergy.XModelica.Media.Water.StandardWater.setState_phX(p,h);

  T:=Exergy.XModelica.Media.Water.StandardWater.temperature(state);
  ex:=Exergy.XModelica.Media.Water.StandardWater.specificExergy(state, refEnv);
  vaperQ:=Exergy.XModelica.Media.Water.StandardWater.vapourQuality(state);

  if vaperQ>0.9999 then
    vaperQ:=1;
  elseif vaperQ<0.0001 then
    vaperQ:=1;
  else
    vaperQ:=vaperQ;
  end if;
  Modelica.Utilities.Streams.print(String(state.p/1e6)+","+String(state.h/1000)+","+String(T-273.15)+","+String(ex/1000)+","+String(vaperQ),fileName);
  //k_mair:=Modelica.Media.Air.MoistAir.k_mair;
  //state.p:=refEnv.p;

 // for T in 273:1:(273+60) loop
 //   for X in 0.0001:0.001:0.042 loop
  //state.T:=T;
  //state.X[1]:=X;
  //state.X[2]:=1-state.X[1];

  //stateCom:=state;

  //ps:=Modelica.Media.Air.MoistAir.saturationPressure(T);
  //omega:=k_mair*ps/(state.p-ps);
  //satX:=omega/(1.0+omega);

 // if X>=satX then
   // stateCom.X[1]:=satX;
   // stateCom.X[2]:=1.0-stateCom.X[1];
  //end if;

   // phi:=Modelica.Media.Air.MoistAir.relativeHumidity(stateCom);
  // ex:=1;
  //phi:=Exergy.XModelica.Media.Air.MoistAir.relativeHumidity_pTX(state.p,state.T,state.X);
  //if phi<0.999 then
 // ex:=Exergy.XBuildings.Media.Air.specificExergy(state, refEnv);
 //outs:=String(T)+","+String(state.X[1])+","+String(ex);

//  Modelica.Utilities.Streams.close("d.txt");

  //else
  //end if;
  //  end for;
  //end for;
end WaterPrint_ph;
