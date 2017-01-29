within Exergy.XModelica.Media.Test;
function specificHumToVaporMolFra
  input Real specificHum;
  output Real vaporMolFra;
  // specificHum is omega
  //    Real molSpecificHum=1.608*specificHum;
  //dryair.MM/steam.MM
protected
   Real molSpecificHum=(1.0/Exergy.XModelica.Media.Air.MoistAir.k_mair)*specificHum;
algorithm

    vaporMolFra:=molSpecificHum/(1 + molSpecificHum);

end specificHumToVaporMolFra;
