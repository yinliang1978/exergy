within Exergy.XClaRa.Components.Mills.HardCoalMills.Fundamentals;
record RollerBowlMillDefinition
  "A record for defining the grey-box model by Niemzcek"

extends ClaRa.Basics.Icons.RecordIcon;
  parameter Real K_1=0.0390 "Fraction of coal on table pulverized per second";
  parameter Real K_2=0.0295
    "Power consumation per kg pulverized coal on table in %.";
  parameter Real K_3=0.0451 "Power consumation per kg raw coal on table in %";
  parameter Real K_4=0.7664
    "Fraction of coal in transport area leaving the mill";
  parameter Real K_5=0.0049
    "Fraction of coal picked ip by primary airper second";
  parameter Real K_6=2.7329 "Classifier speed where m_flow_out=0";
  parameter Real K_7=3.6696
    "Correction factor for primary air difference pressure";
  parameter Real K_8=0.0285 "Nominal value for 'gz/V' of the mill";
  parameter Real K_9=0.5222 "Fraction of coal in transport area falling down";
  parameter Real K_10=5.46 "Fraction of power dissipated";
  parameter Real K_11=19.8e6 "Heat capacity of the mill content";
  parameter Real K_12=1.7
    "Friction loss coefficient for idle load, i.e. zeta/(2A^2)$";
  parameter Real E_e=0.34 "Power consumed for running empy mill in p.u.";
  parameter ClaRa.Basics.Units.Power P_nom=645e3 "Nominal mill power";
end RollerBowlMillDefinition;
