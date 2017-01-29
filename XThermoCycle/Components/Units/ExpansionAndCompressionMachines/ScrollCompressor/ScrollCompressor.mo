within Exergy.XThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor;
model ScrollCompressor "Semi-empirical model of hermetic scroll compressor"
  import ThermoCycle;

parameter Real r_v_in=2.55 "Built-in volume ratio";
parameter Modelica.SIunits.Volume Vs=97.7e-6 "Swept volume";
parameter Modelica.SIunits.ThermalConductance AU_amb=5.382
    "Heat transfer conductance for the ambient heat losses";
parameter Modelica.SIunits.ThermalConductance AU_su=18
    "Heat transfer conductance for the inlet heat exchange";
parameter Modelica.SIunits.ThermalConductance AU_ex=35
    "Heat transfer conductance for the outlet heat exchange";
parameter Modelica.SIunits.Area A_leak=4.5e-7 "Leakage equivalent orifice area";
parameter Modelica.SIunits.Power Wdot_loss_0=242
    "Constant (electro)mechanical losses";
parameter Real alpha=0.2 "Proportionality factor for the proportional losses";
parameter Modelica.SIunits.Length d_ex=0.0075
    "Exhaust pressure drop equivalent orifice diameter";
parameter Modelica.SIunits.Length d_su=0.1
    "Supply pressure drop equivalent orifice diameter";
parameter Modelica.SIunits.Mass m=20 "Total mass of the compressor";
parameter Modelica.SIunits.SpecificHeatCapacity c=466
    "Specific heat capacity of the metal";
parameter Modelica.SIunits.Inertia J=0.02 "Moment of inertia of the rotor";

parameter Modelica.SIunits.Temperature T_su_start = 273.15 + 20
    "Fluid temperature start value, inlet"     annotation (Dialog(tab="Initialization"));
parameter Modelica.SIunits.Temperature T_ex_start = 273.15 + 100
    "Fluid temperature start value, inlet"     annotation (Dialog(tab="Initialization"));
parameter Modelica.SIunits.AbsolutePressure p_su_start = 6e5
    "Fluid pressure start value, inlet"     annotation (Dialog(tab="Initialization"));
parameter Modelica.SIunits.AbsolutePressure p_ex_start = 20e5
    "Fluid pressure start value, outlet"     annotation (Dialog(tab="Initialization"));
parameter Modelica.SIunits.MassFlowRate Mdot_nom=0.13 "Nominal Mass Flow rate"  annotation (Dialog(tab="Initialization"));

  ThermoCycle.Components.Units.PdropAndValves.DP dp_su(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
    constinit=false,
    UseNom=false,
    DELTAp_quad_nom=0,
    use_rho_nom=true,
    UseHomotopy=false,
    A=Modelica.Constants.pi*d_su^2/4,
    p_nom=p_su_start,
    T_nom=T_su_start,
    Mdot_nom=Mdot_nom)
    annotation (Placement(transformation(extent={{-138,-50},{-110,-22}})));
  ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.Isothermal_wall
                  paroi_isotherme(C=m*c, AU_amb=AU_amb)
    annotation (Placement(transformation(extent={{-100,-100},{110,-42}})));
  ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ElectricDrive electricDrive(fstart=50, Np=1)
    annotation (Placement(transformation(extent={{-106,36},{-140,70}})));
  ThermoCycle.Components.Units.HeatExchangers.Semi_isothermal_HeatExchanger
                                HX_su(redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
                                                                             AU=
       AU_su)
    annotation (Placement(transformation(extent={{-100,-52},{-66,-20}})));
  ThermoCycle.Components.Units.HeatExchangers.Semi_isothermal_HeatExchanger
                                HX_ex(redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
                                                                             AU=
       AU_ex)
    annotation (Placement(transformation(extent={{58,-52},{90,-20}})));
  ThermoCycle.Components.Units.PdropAndValves.DP dp_ex(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
    constinit=false,
    UseNom=false,
    UseHomotopy=false,
    A=Modelica.Constants.pi*d_ex^2/4,
    p_nom=p_ex_start,
    T_nom=T_ex_start,
    Mdot_nom=Mdot_nom,
    DELTAp_quad_nom=20000)
    annotation (Placement(transformation(extent={{104,-48},{128,-24}})));
  Nozzle         Internal_leakage(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
    Use_gamma=false,
    A_leak=A_leak,
    P_su_start=p_ex_start,
    T_su_start=T_ex_start,
    P_ex_start=p_su_start)
    annotation (Placement(transformation(extent={{26,-36},{-12,-10}})));
  Volumetric_compression
                compression_4_1(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22,
    r_v_in=r_v_in,
    V_s=Vs,
    T_su_start=T_su_start,
    P_su_start=p_su_start,
    P_ex_start=p_ex_start)
    annotation (Placement(transformation(extent={{-24,-12},{36,40}})));
  Modelica.Mechanics.Rotational.Components.Inertia inertia(
    phi(start=0),
    w(start=300),
    a(start=0),
    J=J)
    annotation (Placement(transformation(extent={{-62,40},{-88,66}})));
  Mechanical_losses mechanical_losses(                         alpha=alpha,
      W_dot_loss_0=Wdot_loss_0)
    annotation (Placement(transformation(extent={{-12,40},{-42,66}})));
  ThermoCycle.Interfaces.Fluid.FlangeA su(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22)
    annotation (Placement(transformation(extent={{-166,-42},{-154,-30}}),
        iconTransformation(extent={{-74,-26},{-66,-18}})));
  ThermoCycle.Interfaces.Fluid.FlangeB ex(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22)
                                            annotation (Placement(transformation(
          extent={{154,-44},{170,-28}}),
                                     iconTransformation(extent={{66,20},{74,28}})));
  Modelica.Blocks.Interfaces.RealInput frequency annotation (Placement(
        transformation(extent={{-172,64},{-132,104}}), iconTransformation(
          extent={{62,-42},{50,-30}})));
  Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a Ambient
    annotation (Placement(transformation(extent={{-4,-104},{12,-88}}),
        iconTransformation(extent={{-4,-90},{4,-82}})));
equation

  connect(compression_4_1.flange_elc, mechanical_losses.flange_A) annotation (
      Line(
      points={{5.7,27.26},{5.7,53},{-11.7,53}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(inertia.flange_a, mechanical_losses.flange_B) annotation (Line(
      points={{-62,53},{-42,53}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(compression_4_1.ex, HX_ex.supply)                   annotation (Line(
      points={{36,14},{46,14},{46,-36},{57.68,-36}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(Internal_leakage.su, HX_ex.supply)                   annotation (Line(
      points={{26.38,-22.74},{46,-22.74},{46,-36},{57.68,-36}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(electricDrive.shaft, inertia.flange_b) annotation (Line(
      points={{-108.38,53},{-88,53}},
      color={0,0,0},
      smooth=Smooth.None));
  connect(dp_su.OutFlow, HX_su.supply)                annotation (Line(
      points={{-111.4,-36},{-100.34,-36}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(paroi_isotherme.W_dot_loss, mechanical_losses.Wdot_loss) annotation (
      Line(
      points={{-58,-65.78},{-58,46.5},{-39.3,46.5}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(paroi_isotherme.Compresseur, HX_ex.port_th)
    annotation (Line(
      points={{3.95,-64.33},{3.95,-58},{74,-58},{74,-51.04}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(paroi_isotherme.Compresseur, HX_su.port_th)
    annotation (Line(
      points={{3.95,-64.33},{3.95,-58},{-83,-58},{-83,-51.04}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(HX_su.exhaust, Internal_leakage.ex)                  annotation (
      Line(
      points={{-66.34,-35.68},{-39.17,-35.68},{-39.17,-22.74},{-12,-22.74}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(HX_su.exhaust, compression_4_1.su)                  annotation (
      Line(
      points={{-66.34,-35.68},{-40,-35.68},{-40,14},{-24,14}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(HX_ex.exhaust, dp_ex.InFlow)                   annotation (Line(
      points={{89.68,-35.68},{96.84,-35.68},{96.84,-36},{105.2,-36}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(dp_ex.OutFlow, ex) annotation (Line(
      points={{126.8,-36},{162,-36}},
      color={0,0,255},
      smooth=Smooth.None));
  connect(frequency, electricDrive.f) annotation (Line(
      points={{-152,84},{-123.68,84},{-123.68,68.98}},
      color={0,0,127},
      smooth=Smooth.None));
  connect(Ambient, paroi_isotherme.Ambient) annotation (Line(
      points={{4,-96},{4,-77.67},{3.95,-77.67}},
      color={191,0,0},
      smooth=Smooth.None));
  connect(su, dp_su.InFlow) annotation (Line(
      points={{-160,-36},{-136.6,-36}},
      color={0,0,255},
      smooth=Smooth.None));
   annotation (Diagram(coordinateSystem(preserveAspectRatio=false, extent={{-160,
            -100},{160,100}}),      graphics={
        Text(
          extent={{-4,-32},{28,-38}},
          lineColor={0,0,0},
          lineThickness=0.5,
          textString="Leak"),
        Text(
          extent={{-100,-16},{-68,-22}},
          lineColor={0,0,0},
          lineThickness=0.5,
          textString="HX1"),
        Text(
          extent={{58,-16},{90,-22}},
          lineColor={0,0,0},
          lineThickness=0.5,
          textString="HX2")}),                 Icon(coordinateSystem(extent={{-160,
            -100},{160,100}}, preserveAspectRatio=false), graphics={
        Ellipse(
          extent={{-60,74},{60,-26}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{-60,30},{60,-90}},
          lineColor={0,0,0},
          fillColor={95,95,95},
          fillPattern=FillPattern.Solid,
          lineThickness=0.5),
        Rectangle(
          extent={{-76,-16},{-60,-30}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{60,30},{76,16}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,255},
          fillPattern=FillPattern.Solid),
        Rectangle(
          extent={{26,-8},{60,-50}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-66,96},{68,74}},
          lineColor={0,0,0},
          lineThickness=0.5,
          fillColor={0,0,0},
          fillPattern=FillPattern.Solid,
          textString="%name")}),
    experiment(StopTime=50),
    __Dymola_experimentSetupOutput,
    Documentation(info="<html>
<p><h4><font color=\"#008000\">Short Description :</font></h4></p>
<p>This model describes a hermetic scroll compressor. The model involves a limited number of parameters physically meaningful. The model accounts for the supply heating-up and exhaust cooling down of the gas, an internal leakage, the internal compression ratio and the electromechanical losses. The model is the one proposed by Winandy et al. (2002). Information used to identify the parameters is extracted from Copeland&apos;s catalogue (ZR72KC-TFD)</p>
<p><br/><h4>Nomenclature :</h4></p>
<p>A                                        Area, m^2</p>
<p>AU                                        Global heat transfer coefficient, W/K</p>
<p>c_p                                        Specific heat, J/kg-K</p>
<p>corr                                        Correction factor for the swept volume, -</p>
<p>d                                        Diameter, m</p>
<p>h                                        Specific enthalpy, J/kg</p>
<p>M_dot                                        Mass flow rate        , kg/s</p>
<p>NTU                                        Number of transfer units,-</p>
<p>N_rot                                        Rotational speed, rpm</p>
<p>p                                        Pressure, Pa</p>
<p>Q_dot                                        Thermal power, W</p>
<p>R                                        Thermal resistance, K/W</p>
<p>r_p                                        Pressure ratio,-</p>
<p>r_v_in                                        Internal built-in volume ratio,-</p>
<p>s                                        Specific entropy, J/kg-K</p>
<p>t                                        Temperature, C</p>
<p>u                                        Specific internal energy, J/kg</p>
<p>v                                        Specific volume, m^3/kg</p>
<p>V_dot_s                                        Swept volume, m^3/s</p>
<p>W_dot                                        Electrical power, W</p>
<p>x                                        quality,-</p>
<p><br/><h4>Greek symbols</h4></p>
<p>alpha                        Proportionality factor for electro-mechanical losses proportional to the internal power</p>
<p>DELTAt                        Temperature difference, K</p>
<p>epsilon                        Efficiency,-</p>
<p>gamma                        Isentropic coefficient,-</p>
<p>m3h                        Cubic meters per hour</p>
<p>rho                        Density, kg m-3</p>
<p><br/><h4>Subscripts</h4></p>
<p>amb                                Ambient</p>
<p>calc                                Calculated</p>
<p>cd                                Condenser</p>
<p>cp                                Compressor</p>
<p>crit                                Critic</p>
<p>ev                                Evaporator</p>
<p>ex                                Exhaust</p>
<p>ex2                                Exhaust before cooling exhaust down</p>
<p>ex1                                Exhaust before exhaust pressure drop</p>
<p>n                                Internal</p>
<p>in*                                Corrected internal</p>
<p>leak                                Leakage</p>
<p>loss                                Electro-mechanical losses</p>
<p>loss0                                Constant electro-mechanical losses</p>
<p>man                                Manufacturer (in Figures)</p>
<p>n                                Nominal</p>
<p>oh                                Over-heating</p>
<p>r                                Refrigerant</p>
<p>s                                Isentropic</p>
<p>s                                Swept volume</p>
<p>sat                                Saturation</p>
<p>sc                                Sub-cooling</p>
<p>su                                Supply</p>
<p>su1                                Supply after supply heating-up</p>
<p>su2                                Supply after mixing with internal mixing</p>
<p>thr                                Throat</p>
<p>w                                Water, enveloppe</p>
<p><br/><h4>References :</h4></p>
<p>[1] Winandy, E., C., Saavedra O., J., Lebrun (2002) Experimental analysis and simplified modelling of a hermetic scrol refrigeration compressor. Applied thermal Engineering 22, 107-120.</p>
<p><br/><h4>Disclaimer</h4></p>
<p>the accuracy or reliability of information presented in this model is not guaranteed or warranted in any way. Every use of this model, for commercial purpose or not, incurs liability of the user only. this model is freely distributed and may not be sold or distributed for commercial purpose. the user is asked to cite his sources and the origin of this model.</p>
<p>!Help us improving this model : any feedback comment is very welcome </p>
<p><br/>Date : January 2014</p>
<p>Authors: Sylvain Quoilin</p>
<p>University of Li&egrave;ge</p>
<p>Faculty of Applied Sciences</p>
<p>thermodynamics Laboratory</p>
<p>Campus of Sart-tilman, B49 (P33)</p>
<p>B-4000 LIEGE (BELGIUM)</p>
<p>Contact: squoilin@ulg.ac.be</p>
<p>website : www.labothap.ulg.ac.be</p>
</html>"));
end ScrollCompressor;
