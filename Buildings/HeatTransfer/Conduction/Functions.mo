within Buildings.HeatTransfer.Conduction;
package Functions
  extends Modelica.Icons.Package;
  function delta_L

    input Modelica.SIunits.Temperature T;
    constant Modelica.SIunits.Pressure p_L=101300;

    output Real value;
  algorithm

    value := 2.0e-7 * (T^0.81) /p_L;

    annotation (preferredView="info", Documentation(info="<html>
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
    value := Modelica.Math.exp(C1/T + C2 + C3*T + C4*T^2 + C5*T^3 + C6*
      Modelica.Math.log(T));
    annotation (preferredView="info", Documentation(info="<html>
  Function to compute the vapour pressure from 273,15K to 473,15K given the absolute temperature (ASHRAE 2001). 
<br/>

</html>"));
  end p_sat;

  function dw_dphi "Function calculating the derivative relative humidity"
    input Buildings.HeatTransfer.Conduction.BaseClasses.RelativeHumidity phi;
    input Modelica.SIunits.MassConcentration w_f;
    input Real b;
    input Real tab_layer[:, :];
    input Real por;
    input Integer switch;
    parameter Integer n=size(tab_layer, 1);
    parameter Real phi_max=1.01;
    Real temp_value[n, 2];
    Real w_max;
    constant Modelica.SIunits.Density Rho=1000;
    output
      Buildings.HeatTransfer.Conduction.BaseClasses.MoistureStorageCapacity
      dw_dphi;

    constant Real B_H=1.0105;
  algorithm
    if (switch == 1) then
      dw_dphi := w_f*(b - 1)*b/((b - phi)^2);

    elseif (switch == 2) then
      temp_value := Buildings.HeatTransfer.Conduction.Functions.der_sorp(
        tab_layer);
      dw_dphi := Modelica.Math.Vectors.interpolate(
          temp_value[:, 1],
          temp_value[:, 2],
          phi);

    else
      w_max := por*Rho;
      dw_dphi := (w_max*B_H*(B_H - phi_max)/phi_max)/(B_H - phi)^2;

    end if;

  end dw_dphi;

  function lambdaHM
    input Modelica.SIunits.Density d;
    input Modelica.SIunits.ThermalConductivity k;
    input Modelica.SIunits.MassConcentration w;
    input Real b;
    input Real tab_layer[:, :];
    input Integer switch;

    output Modelica.SIunits.ThermalConductivity lambdaHM;
  algorithm
    if (switch == 1) then
      lambdaHM := k*(1 + b*w/d);

    else
      lambdaHM := Modelica.Math.Vectors.interpolate(
          tab_layer[:, 1],
          tab_layer[:, 2],
          w);
    end if;
  end lambdaHM;

  function der_sorp
    "Function for the calculation of the derivative of w by phi"
    input Real sorp_tab_layer[:, :];
    parameter Integer n=size(sorp_tab_layer, 1);
    Real temp_value[n, 1];
    output Real value[n, 2];

  algorithm
    for j in 1:n - 1 loop
      temp_value[j, 1] := (sorp_tab_layer[j + 1, 2] - sorp_tab_layer[j, 2])/(
        sorp_tab_layer[j + 1, 1] - sorp_tab_layer[j, 1]);
    end for;
    temp_value[n, 1] := 1E-25;
    value := [sorp_tab_layer[:, 1], temp_value[:, 1]];
  end der_sorp;

  function sorption
    input Buildings.HeatTransfer.Conduction.BaseClasses.RelativeHumidity phi;
    input Modelica.SIunits.MassConcentration w_f;
    input Real b;
    input Real tab_layer[:, :];
    input Real por;
    input Integer switch;
    output Real w;
    constant Modelica.SIunits.Density Rho=1000;
    constant Real B_H=1.0105;
    parameter Real phi_max=1.01;
    Modelica.SIunits.MassConcentration w_max;

  algorithm
    if (switch == 1) then
      w := (w_f*(b - 1)*phi)/(b - phi);
    elseif (switch == 2) then
      w := Modelica.Math.Vectors.interpolate(
          tab_layer[:, 1],
          tab_layer[:, 2],
          phi);
    else
      w_max := por*Rho;
     w := (-w_max*B_H*(B_H - phi_max))/B_H + (w_max*B_H*(B_H - phi_max))/(B_H - phi);

    end if;

  end sorption;

  function Dw
    input Modelica.SIunits.MassConcentration w;
    input Modelica.SIunits.MassConcentration w_f;
    input Real A_layer;
    input Real tab_dww[:, :];
    input Real tab_dws[:, :];
    input Integer switch;
    input Boolean activatesuction;

    output
      Buildings.HeatTransfer.Conduction.BaseClasses.CapillaryTransportCoefficient
      Dw;

  algorithm
    if (switch == 1) then
      if activatesuction then
        Dw := 3.8*((A_layer/w_f)^2)*(1000^((w/w_f) - 1));
      else
        Dw := (3.8*(A_layer/w_f)^2*1000^(w/w_f - 1))*0.1;
      end if;

    else
      if activatesuction then
        Dw := Modelica.Math.Vectors.interpolate(
            tab_dws[:, 1],
            tab_dws[:, 1],
            w);
      else
        Dw := Modelica.Math.Vectors.interpolate(
            tab_dww[:, 1],
            tab_dww[:, 2],
            w);
      end if;
    end if;

  end Dw;

  function w_Kunz

  input Real w_f;
    input Real w_80;
    input Real b;
    input Real phi;
    output Modelica.SIunits.MassConcentration w;
  algorithm
    w := w_f * (b - 1) * phi / (b - phi);

  end w_Kunz;
end Functions;
