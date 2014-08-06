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
</html>",   revisions="<html>
<ul>
<li>
March 17 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

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
    input Buildings.HeatTransfer.Types.RelativeHumidity phi;
    input Modelica.SIunits.MassConcentration w_f;
    input Real b;
    input Real tab[:, :];
    input Real por;
    input Integer switch;
    parameter Integer n=size(tab, 1);
    parameter Real phi_max=1.01;
    Real temp_value[n, 2];
    Real w_max;
    constant Modelica.SIunits.Density Rho=1000;
    output Buildings.HeatTransfer.Types.MoistureStorageCapacity
      dw_dphi;

    constant Real B_H=1.0105;
  algorithm
    if (switch == 1) then
      dw_dphi := w_f*(b - 1)*b/((b - phi)^2);

    elseif (switch == 2) then
      temp_value := Buildings.HeatTransfer.Conduction.Functions.der_tabl(
        tab);
      dw_dphi := Modelica.Math.Vectors.interpolate(
          temp_value[:, 1],
          temp_value[:, 2],
          phi);

    else
      w_max := por*Rho;
      dw_dphi := (w_max*B_H*(B_H - phi_max)/phi_max)/(B_H - phi)^2;

    end if;
  annotation (preferredView="info", Documentation(info="<html>
Function to compute the moisture storage capacity (derivative of w by &phi;).
<p><table>
<tr>
<td> switch = 1</td>
<td> Kunzel's approximation</td>
<td> parameter : w<sub>f</sub>, w<sub>80</sub></td>
</tr>
<tr>
<td>switch = 2</td>
<td>Interpolation of the derivative of the sorption isotherm       </td>
<td> parameter : table of the sorption isotherm</td>
</tr>
<tr>
<td>switch = 3,4,...      </td>
<td>Approximation of the soption isotherm </td>
<td> parameter : porosity of material </td>
</tr>
</table>
<br/>
</html>",   revisions="<html>
<ul>
<li>
May 5 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

          </html>"));
  end dw_dphi;

  function lambdaHM
    input Modelica.SIunits.Density d;
    input Modelica.SIunits.ThermalConductivity k;
    input Modelica.SIunits.MassConcentration w;
    input Real b;
    input Real tab[:, :];
    input Integer switch;

    output Modelica.SIunits.ThermalConductivity lambdaHM;
  algorithm
    if (switch == 1) then
      lambdaHM := k*(1 + b*w/d);

    else
      lambdaHM := Modelica.Math.Vectors.interpolate(
          tab[:, 1],
          tab[:, 2],
          w);
    end if;
  annotation (preferredView="info", Documentation(info="<html>
Function to compute thermal conductivity in moist building material.
<p><table>
<tr>
<td> switch = 1</td>
<td> Kunzel's approximation</td>
<td> parameter : b : thermal conductivity supplement related to the water content in mass percent</td>
</tr>
<tr>
<td>switch = 2, 3...  </td>
<td>Interpolation of table of the moist thermal    conductivity  </td>
<td> parameter : table of the moist thermal conductivity</td>
</tr>
<tr>

</tr>
</table>
<h4>References</h4>


<p>
Cammerer J. and Achtziger J. 1985. <u>Effect of the moisture content on the thermal cconductivity of building materials and insulation products></u> Franhofer Institute.
</p>
<br/>
</html>",   revisions="<html>
<ul>
<li>
Jun 8 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

          </html>"));
  end lambdaHM;

  function der_tabl
    "Function for the calculation of the derivative of w by phi"

  input Real sorp_tab[:, :];
    parameter Integer n=size(sorp_tab, 1);
    Real temp_value[n, 1];
    output Real value[n, 2];

  algorithm
    for j in 1:n - 1 loop
      temp_value[j, 1] := (sorp_tab[j + 1, 2] - sorp_tab[j, 2])/(
        sorp_tab[j + 1, 1] - sorp_tab[j, 1]);
    end for;
    temp_value[n, 1] := 1E-25;
    value := [sorp_tab[:, 1], temp_value[:, 1]];

    annotation (preferredView="info", Documentation(info="<html>
  This function calculates the derivative of a table (in this package the sorption isotherm table) :
  <math><p align=\"center\">  der_tabl [:,1] = tabl [:,1]</p>
  <p align=\"center\">  der_tabl [i,2] = (tabl[i+1,2]-tabl[i,2])/(tabl[i+1,1]-tabl[i,1])</p></math>
<br/>
</html>",   revisions="<html>
<ul>
<li>
Jun 8 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

          </html>"));
  end der_tabl;

  function sorption
    input Buildings.HeatTransfer.Types.RelativeHumidity phi;
    input Modelica.SIunits.MassConcentration w_f;
    input Real b;
    input Real tab[:, :];
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
          tab[:, 1],
          tab[:, 2],
          phi);
    else
      w_max := por*Rho;
     w := (-w_max*B_H*(B_H - phi_max))/B_H + (w_max*B_H*(B_H - phi_max))/(B_H - phi);

    end if;
  annotation (preferredView="info", Documentation(info="<html>
Function to compute water content given the relative humidity.
<p><table>
<tr>
<td> switch = 1</td>
<td> Kunzel's approximation</td>
<td> parameter : w<sub>f</sub>, w<sub>80</sub></td>
</tr>
<tr>
<td>switch = 2</td>
<td>Interpolation of the sorption isotherm       </td>
<td> parameter : table of the sorption isotherm</td>
</tr>
<tr>
<td>switch = 3,4,...      </td>
<td>Approximation of the soption isotherm </td>
<td> parameter : porosity of material </td>
</tr>
</table>
<br/>
</html>",   revisions="<html>
<ul>
<li>
Jun 8 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

          </html>"));
  end sorption;

  function Dw
    input Modelica.SIunits.MassConcentration w;
    input Modelica.SIunits.MassConcentration w_f;
    input Real A_layer;
    input Real tab_dww[:, :];
    input Real tab_dws[:, :];
    input Integer switch;
    input Boolean activatesuction;

    output Buildings.HeatTransfer.Types.CapillaryTransportCoefficient
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
            tab_dws[:, 2],
            w);
      else
        Dw := Modelica.Math.Vectors.interpolate(
            tab_dww[:, 1],
            tab_dww[:, 2],
            w);
      end if;
    end if;

  annotation (preferredView="info", Documentation(info="<html>
This function calculates the liquid transport coefficient for suction D<sub>ws</sub> or for capillary redistribution D<sub>ww</sub>. If the boolean input <var>activatsuction</var> is true 
the material is in suction phase and D<sub>ws</sub> is computed else D<sub>ww</sub> is computed.

<p>
If switch_D<sub>w</sub> = 1
<p align=\"center\" style=\"font-style:italic;\"> D<sub>w</sub> = 3.8&sdot;(A/w<sub>f</sub>)<sup>2</sup>&sdot;1000<sup>w/w<sub>f</sub>-1</p>
Else D<sub>w</sub> is interpolated from the table of the material record.

<h4>References</h4>


<p>
Kiessl, K. 1980. <u>Capillary and vaporous moisture transport in multi-layered building components.</u> Essen University press.
</p>
<br/>
</html>",   revisions="<html>
<ul>
<li>
Jun 14 2014, by Antoine Hudzik:<br/>

First implementation.
</li>
</ul>

          </html>"));
  end Dw;

end Functions;
