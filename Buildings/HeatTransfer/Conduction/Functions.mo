within Buildings.HeatTransfer.Conduction;
package Functions
  extends Modelica.Icons.Package;
  function delta_L

    input Modelica.SIunits.Temperature T;
    constant Modelica.SIunits.Pressure p_L=101300;
    output Real value;
  algorithm
    value := 2.0e-7*T^0.81/p_L;
   annotation (preferredView="info",
    Documentation(info="<html>
  Function to compute the water vapour diffusion coefficient in air given the absolute temperature and the ambiant air pressure 
<br/>

</html>"));
  end delta_L;

  function p_sat
    // Saturation water pressure [Pa]
    input Modelica.SIunits.Temperature T;

    constant Real C1=-5.8E3;
    constant Real C2=1.391;
    constant Real C3=-4.864E-2;
    constant Real C4=4.176E-5;
    constant Real C5=-1.445E-8;
    constant Real C6=6.545;
  output Modelica.SIunits.Pressure value;

  algorithm
    value := Modelica.Math.exp(C1/T+C2+C3*T+C4*T^2+C5*T^3+C6*Modelica.Math.log(T));
     annotation (preferredView="info",
    Documentation(info="<html>
  Function to compute the vapour pressure from 273,15K to 473,15K given the absolute temperature (ASHRAE 2001). 
<br/>

</html>"));
  end p_sat;
end Functions;
