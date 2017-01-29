within Exergy.XModelica.Media.Test;
function specificHumToVaporMassFra
  input Real specificHum;
  output Real vaporMassFra;

algorithm
    vaporMassFra:=specificHum/(1 + specificHum);

end specificHumToVaporMassFra;
