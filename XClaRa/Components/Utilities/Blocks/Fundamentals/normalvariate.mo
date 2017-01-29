within Exergy.XClaRa.Components.Utilities.Blocks.Fundamentals;
function normalvariate "Normally distributed random variable"
  extends ClaRa.Basics.Icons.Function;
  input Real mu "Mean value";
  input Real sigma "Standard deviation";
  input Real si[3] "Input random seed";
  output Real x "Gaussian random variate";
  output Real so[3] "Output random seed";
protected
  constant Real NV_MAGICCONST=4*exp(-0.5)/sqrt(2.0);
  Real s1[3];
  Real  s2[3];
  Real z;
  Real zz;
  Real u1;
  Real u2;
  Boolean break_=false;
algorithm
  s1 := si;
  u2 := 1;
  while not break_ loop
    (u1,s2) := random(s1);
    (u2,s1) := random(s2);
    z := NV_MAGICCONST*(u1-0.5)/u2;
    zz := z*z/4.0;
    break_ := zz <= (- Math.log(u2));
  end while;
    x := mu + z*sigma;
    so := s1;
end normalvariate;
