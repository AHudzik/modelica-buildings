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

  function dw_dphi "Function calculating the derivative relative humidity"
    input Real      phi;
    input Real      w_f;
    input Real      w_80;
    input Real      b;
    input Real tab_layer[:,:];
    input Boolean Kunzel;
    parameter Integer n = size(tab_layer,1);
    Real temp_value[  n,2];
    output Real dw_dphi;

  algorithm
    if Kunzel==true then
    dw_dphi :=w_f*(b - 1)*b/((b - phi)^2);

    else
      temp_value :=Buildings.HeatTransfer.Conduction.Functions.der_sorp(
        tab_layer);
      dw_dphi:=Modelica.Math.Vectors.interpolate(temp_value[1,:], temp_value[2, :], phi);
    end if;

  end dw_dphi;

  function w_Kunzel
    input Real      w_f;
    input Real      w_80;
    input Real      b;
    input Real      phi;

    output Modelica.SIunits.MassConcentration w;

  algorithm
     w:= w_f*((b - 1)*phi)/(b - phi);

  end w_Kunzel;

  function lambdaHM
    input Modelica.SIunits.Density d;
    input Modelica.SIunits.ThermalConductivity k;
    input Modelica.SIunits.MassConcentration w;
    input Real b;
    input Real tab_layer[ :,:];
    input Boolean Kunzel;

    output Modelica.SIunits.ThermalConductivity lambdaHM;
  algorithm
    if Kunzel == true then
    lambdaHM :=k*(1 + b*w/d);

    else
    lambdaHM:= Modelica.Math.Vectors.interpolate(
        tab_layer[1, :],
        tab_layer[2, :],
        w);
    end if;
  end lambdaHM;

  function der_sorp
    "Function for the calculation of the derivative of w by phi"
  input Real sorp_tab_layer[  :,:];
  parameter Integer n = size(sorp_tab_layer,1);
      Real temp_value[n,1];
      output Real value[  n,2];

  algorithm
          for j in 1:n-1 loop
            temp_value[j,1] := (sorp_tab_layer[j+1,2]-sorp_tab_layer[j,2])/(sorp_tab_layer[j+1,1]-sorp_tab_layer[j,1]);
          end for;
          temp_value[n, 1] :=1E-25;
          value := [sorp_tab_layer[:,1], temp_value[:,1]];
  end der_sorp;

  function sorption
    input Real phi;
    input Real w_f;
    input Real w_80;
    input Real b;
    input Real tab_layer[ :,:];
    input Boolean Kunzel;
    output Real w;

  algorithm
    if Kunzel== true then
      w:= w_f*((b - 1)*phi)/(b - phi);

    else
   w:=Modelica.Math.Vectors.interpolate(
        tab_layer[2, :],
        tab_layer[1, :],
        phi);

    end if;
  end sorption;

  function Capillary_transp_coef
    input Real w;
    input Real w_f;
    input Modelica.SIunits.Area A;
    input Real tab_layer[ :,:];
    input Boolean Kunzel;

  output Real Dw;

  algorithm
    if Kunzel== true then
      Dw:= 3.8*(A/w_f)^2*1000^(w/(w_f - 1));

    else
      Dw:= Modelica.Math.Vectors.interpolate(
        tab_layer[1, :],
        tab_layer[2, :],
        w);
    end if;

  end Capillary_transp_coef;
end Functions;
