within Exergy.Contours;
model WaterContour_dT
 // Real x;
  Exergy.Utilities.RefEnv refEnv;
  parameter String fileName="WaterData.txt";
  parameter Real[:] vrange=-2:0.1:3;
  parameter Real[:] verange=10 .^vrange;
  parameter Real[:] Trange=-3:0.1:3;
  parameter Real[:] Terange=(10 .^Trange).+273.16;
initial algorithm
 // test.pt(refEnv);
equation

  when terminal() then
    Modelica.Utilities.Files.removeFile(fileName);
    //if the code is put in the function ,the memory will be out
    for v in verange loop
     for T in Terange loop
        WaterPrint_dT(
              refEnv,
              v,
              T,
              fileName);
     end for;
       end for;
  end when;
  //x=1;

end WaterContour_dT;
