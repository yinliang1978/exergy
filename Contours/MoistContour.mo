within Exergy.Contours;
model MoistContour
  /*
  
  // for standard
  Exergy.Utilities.RefEnv refEnv;

 // Exergy.Utilities.RefEnv refEnv;
  parameter String fileName="MoistData.txt";
initial algorithm 
 // test.pt(refEnv);
equation 
  when terminal() then
    Modelica.Utilities.Files.removeFile(fileName);
    //if the code is put in the function ,the memory will be out
  //  for T in (273):1:(273+60) loop
  //   for X in 0.0001:0.001:0.042 loop
       for T in (273.15):1:(273.15+60) loop
     for X in 0.0001:0.0005:0.042 loop
        MoistPrint(
          refEnv,
          T,
          X,
          fileName);
     end for;
       end for;
  end when;
  
  */
/*
   // for beijing summer dT 33.5 WT 26.4 omega 0.019
   Real T=273.15+33.5;
   //Real MoistMassX[2]={0.019/(1+0.019),1/(1+0.019)};
   //Real omega=MoistMassX[1]/MoistMassX[2];
   Real omega=0.019;
   Real vaporMolFra=Exergy.XModelica.Media.Test.specificHumToVaporMolFra(omega);

   Exergy.Utilities.RefEnv refEnv(T=T,vaporMolFra=vaporMolFra);

 // Exergy.Utilities.RefEnv refEnv;
  parameter String fileName="MoistData.txt";
initial algorithm 
 // test.pt(refEnv);
equation 
  when terminal() then
    Modelica.Utilities.Files.removeFile(fileName);
    //if the code is put in the function ,the memory will be out
  //  for T in (273):1:(273+60) loop
  //   for X in 0.0001:0.001:0.042 loop
       for T in (273.15):1:(273.15+60) loop
     for X in 0.0001:0.0005:0.042 loop
        MoistPrint(
          refEnv,
          T,
          X,
          fileName);
     end for;
       end for;
  end when;
  */

     // for beijing   winter dT -9.9 phi 44%
   Real T=273.15-9.9;
   //Real MoistMassX[2]={0.019/(1+0.019),1/(1+0.019)};
   //Real omega=MoistMassX[1]/MoistMassX[2];
   Real omega=0.000718;
   Real vaporMolFra=Exergy.XModelica.Media.Test.specificHumToVaporMolFra(omega);

   Exergy.Utilities.RefEnv refEnv(T=T,vaporMolFra=vaporMolFra);

 // Exergy.Utilities.RefEnv refEnv;
  parameter String fileName="MoistData.txt";
initial algorithm
 // test.pt(refEnv);
equation
  when terminal() then
    Modelica.Utilities.Files.removeFile(fileName);
    //if the code is put in the function ,the memory will be out
  //  for T in (273):1:(273+60) loop
  //   for X in 0.0001:0.001:0.042 loop
       for T in (273.15-10):0.5:(273.15+60) loop
     for X in 0.0001:0.0002:0.042 loop
        MoistPrint(
          refEnv,
          T,
          X,
          fileName);
     end for;
       end for;
  end when;

  /*
  
  
 // Real x;
 //beijing summer dT 33.5 WT 26.4 omega 0.019
 // Real T=273.15+33.5;
 // Real MoistMassX[2]={0.019/(1+0.019),1/(1+0.019)};
 // Exergy.Utilities.RefEnv refEnv(T=T,MoistMassX=MoistMassX);

 //beijing winter dT -9.9 phi 44%
    Modelica.SIunits.Temperature T=273.15-9.9;
   Real omega=0.000718;
  Real MoistMassX[2]={omega/(1+omega),1/(1+omega)};
  Exergy.Utilities.RefEnv refEnv(T=T,MoistMassX=MoistMassX);

 // Exergy.Utilities.RefEnv refEnv;
  parameter String fileName="MoistData.txt";
initial algorithm 
 // test.pt(refEnv);
equation 
  when terminal() then
    Modelica.Utilities.Files.removeFile(fileName);
    //if the code is put in the function ,the memory will be out
  //  for T in (273):1:(273+60) loop
  //   for X in 0.0001:0.001:0.042 loop
       for T in (273.15-9.9):1:(273+50) loop
     for X in 0.00001:0.0001:0.042 loop
        MoistPrint(
          refEnv,
          T,
          X,
          fileName);
     end for;
       end for;
  end when;
  //x=1;
*/
end MoistContour;
