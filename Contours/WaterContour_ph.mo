within Exergy.Contours;
model WaterContour_ph
 // Real x;
  Exergy.Utilities.RefEnv refEnv;
  parameter String fileName="WaterData.txt";
  parameter Real[:] prange=3:0.02:8;
  parameter Real[:] perange=10 .^prange;
  parameter Real[:] hrange=5:0.01:6.56;
  parameter Real[:] herange=(10 .^hrange);
initial algorithm
 // test.pt(refEnv);
equation

  when terminal() then
    Modelica.Utilities.Files.removeFile(fileName);
    //if the code is put in the function ,the memory will be out
    for p in perange loop
     for h in herange loop
        WaterPrint_ph(
          refEnv,
          p,
          h,
          fileName);
     end for;
       end for;
  end when;
  //x=1;

end WaterContour_ph;
