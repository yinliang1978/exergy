within Exergy.XThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor;
model Volumetric_compression
  "Model of a volumetric compression including swept volume and under and over-expansion losses"
  import ThermoCycle;
replaceable package Medium = ThermoCycle.Media.DummyFluid;
Medium.ThermodynamicState inlet(s(start=Medium.specificEntropy_pT(P_su_start,T_su_start)));
Medium.ThermodynamicState inside;
Medium.ThermodynamicState outlet;
Medium.Density rho_su(start=Medium.density_pT(P_su_start,T_su_start));
Medium.SpecificEnthalpy h_s;
parameter Medium.Temperature T_su_start = 273.15 + 20
    "Fluid temperature start value, inlet"     annotation (Dialog(tab="Initialization"));
parameter Medium.AbsolutePressure P_su_start = 6e5
    "Fluid pressure start value, inlet"     annotation (Dialog(tab="Initialization"));
parameter Medium.AbsolutePressure P_ex_start = 20e5
    "Fluid pressure start value, outlet"     annotation (Dialog(tab="Initialization"));
Medium.AbsolutePressure P_ad(start = P_su_start*r_v_in,stateSelect = StateSelect.prefer);
parameter Real r_v_in = 2.55;
parameter Modelica.SIunits.Volume V_s=100e-6 "Compressor swept Volume";

Modelica.SIunits.AngularFrequency N_rot "Compressor rotational speed in Hz";
Modelica.SIunits.SpecificEnergy w_1;
Modelica.SIunits.SpecificEnergy w_2;
Modelica.SIunits.Power W_dot_in;
Modelica.SIunits.AngularVelocity omega_m
    "Angular velocity of the shaft [rad/s] ";

  ThermoCycle.Interfaces.Fluid.FlangeA su(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22)
    annotation (Placement(transformation(extent={{-110,-10},{-90,10}})));
  ThermoCycle.Interfaces.Fluid.FlangeB ex(
    redeclare package Medium =
        ThermoCycle.Components.Units.ExpansionAndCompressionMachines.ScrollCompressor.R22)
                                            annotation (Placement(transformation(
          extent={{90,-10},{110,10}}),
                                     iconTransformation(extent={{90,-10},{110,10}})));
  Modelica.Mechanics.Rotational.Interfaces.Flange_b flange_elc annotation (
      Placement(transformation(extent={{-20,32},{12,64}}),iconTransformation(
          extent={{-10,42},{8,60}})));
equation
//Flow rate calculation
su.m_flow = N_rot * V_s * rho_su;

//Thermodynamic states
inlet = Medium.setState_ph(su.p,inStream(su.h_outflow));
outlet = Medium.setState_ph(ex.p,ex.h_outflow);

//Using preferably functions for which the derivative is defined:
rho_su = Medium.density_ph(su.p,inStream(su.h_outflow));

//Internal thermodynamic state
inside = Medium.setState_ps(P_ad,inlet.s);

//After the isentropic expansion
h_s = Medium.specificEnthalpy_ps(P_ad,inlet.s);
rho_su*r_v_in = Medium.density_ph(P_ad,h_s);

//Flow rate:
su.m_flow + ex.m_flow = 0;

//Enthalpies
su.h_outflow = inStream(ex.h_outflow); //by pass in case of flow reversal
ex.h_outflow = inlet.h + w_1 + w_2;

//Isentropic compression
w_1 = (h_s-inlet.h);

//Constant-volume compression
w_2 = (ex.p - P_ad) / inside.d;

//Internal compression power:
W_dot_in = su.m_flow*(w_1 + w_2);

//Mechanical port
omega_m = der(flange_elc.phi);
omega_m = 2*N_rot*Modelica.Constants.pi;
flange_elc.tau = W_dot_in/(2*N_rot*Modelica.Constants.pi);

  annotation (Icon(coordinateSystem(preserveAspectRatio=false, extent={{-100,-100},
            {100,100}}), graphics={
        Polygon(
          points={{-100,-80},{20,-40},{100,-40},{100,40},{20,40},{-100,80},{-100,
              -80}},
          lineColor={0,0,255},
          smooth=Smooth.None,
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid),
        Text(
          extent={{-27,17},{27,-17}},
          lineColor={0,0,255},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          textString="s=Cst",
          origin={-43,-3},
          rotation=0),
        Text(
          extent={{-27,17},{27,-17}},
          lineColor={0,0,255},
          fillColor={170,213,255},
          fillPattern=FillPattern.Solid,
          origin={51,-3},
          rotation=0,
          textString="V=Cst")}));
end Volumetric_compression;
