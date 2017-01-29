within Exergy.XClaRa.Components.Furnace.ChemicalReactions;
model CoalReactionZone
//___________________________________________________________________________//
// Component of the ClaRa library, version: 1.1.2                        //
//                                                                           //
// Licensed by the DYNCAP/DYNSTART research team under Modelica License 2.   //
// Copyright © 2013-2016, DYNCAP/DYNSTART research team.                     //
//___________________________________________________________________________//
// DYNCAP and DYNSTART are research projects supported by the German Federal //
// Ministry of Economic Affairs and Energy (FKZ 03ET2009/FKZ 03ET7060).      //
// The research team consists of the following project partners:             //
// Institute of Energy Systems (Hamburg University of Technology),           //
// Institute of Thermo-Fluid Dynamics (Hamburg University of Technology),    //
// TLK-Thermo GmbH (Braunschweig, Germany),                                  //
// XRG Simulation GmbH (Hamburg, Germany).                                   //
//___________________________________________________________________________//

extends Exergy.XClaRa.Components.Furnace.ChemicalReactions.PartialReactionZone;
//Please note that the following calculations are carried out for 1 kg of fuel! Because the factor of 1 would have no influence it is not shown inside these equations.

equation
  n_C = xi_fuel_in[1]/ClaRa.Basics.Constants.M_C;
  n_H = xi_fuel_in[2]/ClaRa.Basics.Constants.M_H;
  n_O = xi_fuel_in[3]/ClaRa.Basics.Constants.M_O;
  n_N = xi_fuel_in[4]/ClaRa.Basics.Constants.M_N;
                                           // N not N2!!!!!!
  n_S = xi_fuel_in[5]/ClaRa.Basics.Constants.M_S;
  n_Ash = xi_fuel_in[6]/ClaRa.Basics.Constants.M_Ash;
  n_H2O = (1 - sum(xi_fuel_in))/ClaRa.Basics.Constants.M_H2O;

  0 = xi_fuel_in[1]/ClaRa.Basics.Constants.M_C*xi_CO - n_CO;
  0 = xi_fuel_in[4]/ClaRa.Basics.Constants.M_N*xi_NOx - n_NO;
  0 = xi_fuel_in[1]/ClaRa.Basics.Constants.M_C - n_CO2 - n_CO;
  0 = xi_fuel_in[2]/ClaRa.Basics.Constants.M_H/2 + (1 - sum(
    xi_fuel_in))/(ClaRa.Basics.Constants.M_H*2 + ClaRa.Basics.Constants.M_O)
     - n_H2O_prod;
  0 = xi_fuel_in[5]/ClaRa.Basics.Constants.M_S - n_SO2;
  0 = xi_fuel_in[4]/ClaRa.Basics.Constants.M_N/2 - n_N2 - n_NO/2;

 for i in 1:(flueGas.nc-1) loop
     if i==1 then prod_comp[1] =xi_fuel_in[6]*(1 - xi_slag);
     else if i==2 then
        prod_comp[2] = n_CO*(ClaRa.Basics.Constants.M_C + ClaRa.Basics.Constants.M_O);
     else if i==3 then
          prod_comp[3] = n_CO2*(ClaRa.Basics.Constants.M_C + 2*ClaRa.Basics.Constants.M_O);
     else if i==4 then
            prod_comp[4] = n_SO2*(ClaRa.Basics.Constants.M_S + 2*
              ClaRa.Basics.Constants.M_O);
     else if i==5 then
              prod_comp[5] = n_N2*(2*ClaRa.Basics.Constants.M_N);
     else if i==6 then
                prod_comp[6] = -(n_CO/2 + n_CO2 + n_NO/2 + n_H/4.0 +
                  n_S - n_O/2)*ClaRa.Basics.Constants.M_O*2.0;
     else if i==7 then
                  prod_comp[7] = n_NO*(ClaRa.Basics.Constants.M_N +
                    ClaRa.Basics.Constants.M_O);
     else if i==8 then
                    prod_comp[8] = n_H2O_prod*(ClaRa.Basics.Constants.M_H
                      *2 + ClaRa.Basics.Constants.M_O);
    else
      prod_comp[i] = 0;
     end if; end if; end if; end if; end if; end if; end if; end if;
end for;

end CoalReactionZone;
